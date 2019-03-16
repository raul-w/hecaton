#!/usr/bin/env python3

"""
Postprocess VCF files, so that dispersed duplication variants are represented
properly
"""

import argparse
import pysam
import sys
from collections import defaultdict
from itertools import count
from classes.breakend import Breakend
from classes.genomic_coordinate import GenomicCoordinate
from classes.sv import SV

def string2bool(v):
    """
    Change string to Boolean

    See https://stackoverflow.com/questions/15008758/parsing-boolean-values-with-argparse

    :param v: A string that represents a Boolean value
    :return: Boolean
    """
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def parse_cl_args(in_args):
    """
    Parse command line arguments

    :param in_args: All command line arguments
    :return: None
    """
    description = "Parse VCF files, extracting simple CNVs in BEDPE format"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-i", "--input_vcf", type=str,
                        help="Path to VCF file")
    parser.add_argument("-o", "--output_bedpe", type=str,
                        help="Name of BEDPE output file")
    parser.add_argument("-t", "--tool", type=str,
                        help="Tool used to call variants")
    parser.add_argument("-d", "--dispersed_duplications", type=string2bool, default=True,
                        help="True if you want to check for dispersed duplicatons, otherwise False")
    args = parser.parse_args(in_args)
    return args

def parse_breakends_and_insertions_from_vcf(vcf: pysam.VariantFile):
    """
    Parse breakends from vcf

    :param record: A VariantFile instance
    :return: List of Breakend instances
    """
    breakends = {}
    # breakends: key id, value: breakend
    id_generator = count(1)
    # keep track of adjacencies and BNDs associated with them, so we can
    # find mates of BND entries
    adjacencies_dict = defaultdict(int)
    # adjacencies_dict: key (chromA, posA, chromB, posB), value id
    for record in vcf.fetch():
        # initialize inserted sequence variable
        inserted_sequence = "."
        # get breakends according to type
        variant_type = record.info["SVTYPE"]
        if variant_type == "BND":
            # extract mate coordinate from ALT field
            alts = record.alts
            # TODO: allow BNDs to have multiple mates
            alt_field = alts[0]
            # check if [ or ] is in alt field, use it to determine mate direction
            if "[" in alt_field:
                alt_contents = alt_field.split("[")
                mate_direction = "R"
            elif "]" in alt_field:
                alt_contents = alt_field.split("]")
                mate_direction = "L"
            else:
                # breakend with unknown mate, skip it
                continue
            # alt_contents is either [REF, CHROM:POS , ""] or ["", CHROM:POS , REF]
            # ref has a length longer than 1 in case of an insertion
            # determine position of mate breakend
            if alt_contents[0]:
                alt = alt_contents[0]
                # take insertions into account
                if len(alt) > 1:
                    inserted_sequence = alt[1:]
                mate_position = "R"
            else:
                alt = alt_contents[2]
                # take insertions into account
                if len(alt) > 1:
                    inserted_sequence = alt[:-1]
                mate_position = "L"
            # extract genomic position of breakend and mate
            breakend_chrom = record.chrom
            breakend_pos = record.pos
            mate_chrom = alt_contents[1].split(":")[0]
            mate_pos = int(alt_contents[1].split(":")[1])
            breakend_coordinate = GenomicCoordinate(breakend_chrom, breakend_pos)
            mate_coordinate = GenomicCoordinate(mate_chrom, mate_pos)
            new_id = next(id_generator)
            # create breakend
            breakend = Breakend(new_id, breakend_coordinate, mate_coordinate, mate_direction,
                                mate_position, record, inserted_sequence=inserted_sequence)
            # add adjacency to adjacency dict
            breakend_adjacency = (breakend_chrom, breakend_pos, mate_chrom, mate_pos)
            adjacencies_dict[breakend_adjacency] = new_id
            # check if mate is already in breakend dictionary
            mate_adjacency = (mate_chrom, mate_pos, breakend_chrom, breakend_pos)
            mate_id = adjacencies_dict[mate_adjacency]
            breakends[new_id] = breakend
            if mate_id != 0:
                # update mate ids
                breakend.mate_id = mate_id
                breakends[mate_id].mate_id = new_id
        elif variant_type == "DEL":
            left_coordinate = GenomicCoordinate(record.chrom, record.pos)
            right_coordinate = GenomicCoordinate(record.chrom, record.stop + 1)
            left_id = next(id_generator)
            right_id = next(id_generator)
            # add two breakends for each coordinate
            left_breakend = Breakend(left_id, left_coordinate,
                                     right_coordinate, "R", "R", record,
                                     mate_id=right_id)
            right_breakend = Breakend(right_id, right_coordinate,
                                      left_coordinate, "L", "L", record,
                                      mate_id=left_id)
            # add ids
            breakends[left_id] = left_breakend
            breakends[right_id] = right_breakend
        elif variant_type == "DUP":
            left_coordinate = GenomicCoordinate(record.chrom, record.pos + 1)
            right_coordinate = GenomicCoordinate(record.chrom, record.stop)
            left_id = next(id_generator)
            right_id = next(id_generator)
            # add two breakends for each coordinate
            left_breakend = Breakend(left_id, left_coordinate, right_coordinate, "L", "L",
                                     record, mate_id=right_id)
            right_breakend = Breakend(right_id, right_coordinate, left_coordinate, "R", "R",
                                      record, mate_id=left_id)
            breakends[left_id] = left_breakend
            breakends[right_id] = right_breakend
        elif variant_type == "INS":
            inserted_sequence = "N"
            # get inserted sequence, if present
            alts = record.alts
            alt_field = alts[0]
            if alt_field != "<INS>":
                    inserted_sequence = alt_field[1:]
            left_coordinate = GenomicCoordinate(record.chrom, record.pos)
            right_coordinate = GenomicCoordinate(record.chrom, record.pos + 1)
            left_id = next(id_generator)
            right_id = next(id_generator)
            # add two breakends for each coordinate
            left_breakend = Breakend(left_id, left_coordinate,
                                     right_coordinate, "R", "R", record,
                                     mate_id=right_id,
                                     inserted_sequence=inserted_sequence)
            right_breakend = Breakend(right_id, right_coordinate,
                                      left_coordinate, "L", "L", record,
                                      mate_id=left_id,
                                      inserted_sequence=inserted_sequence)
            breakends[left_id] = left_breakend
            breakends[right_id] = right_breakend
        else:
            # not a copy number variant, continue
            continue
    return breakends

def cluster_breakends(breakends: dict):
    """
    Cluster breakends in a dictionary of breakends

    :param breakends: Dictionary with a tuples (chr, pos) as keys and breakends
    as values
    :return: Dictionary of sets of breakends
    """
    cluster_id_generator = count(1)
    # use dictionarys to keep track of cluster_ids of breakends
    breakend_id_dict = defaultdict(int)
    # give empty dictionary if there are no SVs
    if not breakends:
        return breakend_id_dict
    # generate list of breakends, sorted by chromosome and position
    sorted_breakends = list(breakends.values())
    sorted_breakends.sort(
        key=lambda x: [x.coordinate.chromosome, x.coordinate.position])
    # breakend_id_dict: breakend_id as key, cluster_id as value
    cluster_id_dict = defaultdict(set)
    # cluster_dict: cluster_id as key, set of breakpoint_ids
    # case i = 0
    cluster_id = next(cluster_id_generator)
    breakend_id = sorted_breakends[0].identifier
    breakend_mate_id = sorted_breakends[0].mate_id
    breakend_id_dict[breakend_id] = cluster_id
    breakend_id_dict[breakend_mate_id] = cluster_id
    cluster_id_dict[cluster_id].add(breakend_id)
    cluster_id_dict[cluster_id].add(breakend_mate_id)
    # loop through rest of breakends
    for i in range(1, len(sorted_breakends)):
        breakend = sorted_breakends[i]
        # check if breakend is already clustered
        breakend_id = breakend.identifier
        current_cluster_id = breakend_id_dict[breakend_id]
        if current_cluster_id == 0:
            # create new cluster for itself and its mate
            new_cluster_id = next(cluster_id_generator)
            breakend_id_dict[breakend_id] = new_cluster_id
            breakend_mate_id = breakend.mate_id
            breakend_id_dict[breakend_mate_id] = new_cluster_id
            # update cluster id dict
            cluster_id_dict[new_cluster_id].add(breakend_id)
            cluster_id_dict[new_cluster_id].add(breakend_mate_id)
            # update cluster id
            current_cluster_id = new_cluster_id
        # if breakend is within 10 bp downstream of other breakend,
        # put it in same cluster together with its mate
        previous_breakend = sorted_breakends[i - 1]
        previous_breakend_id = previous_breakend.identifier
        if breakend.coordinate.chromosome == previous_breakend.coordinate.chromosome \
                and breakend.coordinate.position - \
                previous_breakend.coordinate.position <= 10:
            new_cluster_id = breakend_id_dict[previous_breakend_id]
            # merge clusters in cluster id dict
            cluster_id_dict[new_cluster_id] = \
                cluster_id_dict[new_cluster_id].union(
                    cluster_id_dict[current_cluster_id])
            # update ids of breakends in breakend dict
            for breakend_id in cluster_id_dict[current_cluster_id]:
                breakend_id_dict[breakend_id] = new_cluster_id
            # delete old cluster id from cluster_id dict
            del cluster_id_dict[current_cluster_id]
        else:
            continue
    # convert the dictionary of sets of breakend coordinates
    # in a dictionary of cluster ids as keys and a set of breakends as values
    breakend_dict = defaultdict(set)
    for breakend_id in breakends:
        breakend = breakends[breakend_id]
        cluster_id = breakend_id_dict[breakend_id]
        breakend_dict[cluster_id].add(breakend)
    return breakend_dict

def positive_strand(sv):
    """
    Change strandA of the structural variant to +, if not already

    :param sv: SV object
    :return: SV object in which strandA and strandB are both +
    """
    if sv.strandA == "-":
        if sv.strandB == "+":
            strandB = "-"
        else:
            strandB = "+"
        new_sv = SV(
            chromosomeA=sv.chromosomeB,
            startA=sv.startB,
            endA=sv.endB,
            chromosomeB=sv.chromosomeA,
            startB=sv.startA,
            endB=sv.endA,
            strandA="+", strandB=strandB,
            sv_type="BND", inserted_sequence=sv.inserted_sequence)
    else:
        new_sv = sv
    return new_sv

def get_bnd_from_breakend(breakend):
    """
    Obtain BND structural variant from breakend object

    :param breakend: Breakend object
    :return: Instance of SV classes of type BND
    """
    # get orientation of strand A from position of mate
    if breakend.mate_position == "R":
        strand_A = "+"
    else:
        strand_A = "-"
    # get orientation of strand B from direction of mate
    if breakend.mate_direction == "R":
        strand_B = "+"
    else:
        strand_B = "-"
    structural_variant = SV(chromosomeA=breakend.coordinate.chromosome,
                              startA=breakend.coordinate.position,
                              endA=breakend.coordinate.position,
                              chromosomeB=breakend.mate_coordinate.chromosome,
                              startB=breakend.mate_coordinate.position,
                              endB=breakend.mate_coordinate.position,
                              strandA=strand_A, strandB=strand_B,
                              sv_type="BND",
                            inserted_sequence=breakend.inserted_sequence)
    return structural_variant

def get_simple_sv_from_breakend(breakend):
    """
    Obtain simple SV (DEL, DUP:TANDEM or INS) from breakend object if possible.

    :param breakend: Breakend object
    :return: Instance of SV classes of type DEL, DUP:TANDEM, INS, or BND
    """
    structural_variant = get_bnd_from_breakend(breakend)
    if structural_variant.chromosomeA != structural_variant.chromosomeB or \
            structural_variant.strandA != structural_variant.strandB:
        # Type of SV with 2 breakends that is not a deletion or tandem duplication
        return structural_variant
    # put the structural variant in + orientation, if not already
    structural_variant = positive_strand(structural_variant)
    chrom = structural_variant.chromosomeA
    if abs(structural_variant.startA - structural_variant.startB) == 1:
        # insertion
        start_pos = structural_variant.startA
        end_pos = structural_variant.startB
        # correct GRIDSS misinterpretation of insertion positions
        if end_pos == start_pos - 1:
            end_pos = start_pos + 1
        inserted_sequence = structural_variant.inserted_sequence
        structural_variant = SV(chromosomeA=chrom,
                                startA=start_pos, endA=end_pos,
                                sv_type="INS",
                                inserted_sequence=inserted_sequence)
    elif structural_variant.startA < structural_variant.startB:
        # deletion
        start_pos = structural_variant.startA + 1
        end_pos = structural_variant.startB - 1
        structural_variant = SV(chromosomeA=chrom,
                                startA=start_pos, endA=end_pos,
                                sv_type="DEL")
    else:
        # tandem duplication
        start_pos = structural_variant.startB
        end_pos = structural_variant.startA
        structural_variant = SV(chromosomeA=chrom,
                                startA=start_pos, endA=end_pos,
                                sv_type="DUP:TANDEM")
    return structural_variant

def add_features(sv, variant_record):
    """
    Add features from a variant record to the sv

    :param sv: SV object
    :param variant_record: VariantRecord object from pysam
    :return: SV object containing features
    """
    # LUMPY format
    if "PE" in variant_record.format:
        sv.read_pairs_defined = 1
        sv.read_pairs = str(variant_record.samples[0]["PE"])
    # Manta format
    elif "PR" in variant_record.format:
        sv.manta_ref_read_pairs_defined = 1
        sv.read_pairs_defined = 1
        sv.manta_ref_read_pairs = \
        variant_record.samples[0]["PR"][0]
        sv.read_pairs = str(variant_record.samples[0]["PR"][1])
    # GRIDSS format
    elif "RP" in variant_record.format:
        sv.read_pairs_defined = 1
        sv.read_pairs = str(variant_record.samples[0]["RP"])
    # Delly format
    elif "PE" in variant_record.info:
        sv.read_pairs_defined = 1
        sv.read_pairs = str(variant_record.info["PE"])
    # LUMPY, Manta and GRIDSS format
    if "SR" in variant_record.format:
        # Manta format
        if isinstance(variant_record.samples[0]["SR"], tuple):
            sv.manta_ref_split_reads_defined = 1
            sv.manta_ref_split_reads = \
                variant_record.samples[0]["SR"][0]
            sv.split_reads_defined = 1
            sv.split_reads = str(variant_record.samples[0]["SR"][1])
        else:
            sv.split_reads_defined = 1
            sv.split_reads = str(variant_record.samples[0]["SR"])
    # DELLY format
    elif "SR" in variant_record.info:
        sv.split_reads_defined = 1
        sv.split_reads = str(variant_record.info["SR"])
    if "LowQual" in variant_record.filter:
        sv.delly_low_qual = 1
    return sv

def structural_variant_from_breakend_cluster(breakend_cluster, check_dispersed_duplication=True):
    """
    Extract structural variant from break-end cluster

    :param breakend_cluster: Set of breakends
    :param check_dispersed_duplication: Will check for dispersed duplications
    if True (boolean)
    :return: List of structural variants
    """
    structural_variants = []
    variant_records = []
    mate_ids = []
    # sort breakends by chromosome and position
    breakend_cluster = list(breakend_cluster)
    breakend_cluster.sort(key=lambda x: [x.coordinate.chromosome, x.coordinate.position])
    if len(breakend_cluster) == 4:
        # get the first breakend
        breakend_a = breakend_cluster[0]
        structural_variant_a = get_bnd_from_breakend(breakend_a)
        mate_ids.append(breakend_a.mate_id)
        # get the variant record
        variant_record_a = breakend_a.variant_record
        # get the second breakend
        breakend_b = breakend_cluster[1]
        if breakend_b.identifier in mate_ids:
            # breakend is mate of first breakend, go to next one
            breakend_b = breakend_cluster[2]
        structural_variant_b = get_bnd_from_breakend(breakend_b)
        # get the variant record
        variant_record_b = breakend_b.variant_record
        # check if the two variants add up to a dispersed duplication
        # if strand A != strand B in any of the two, we have no dispersed duplications,
        # but a structural variant involving an inversion
        if structural_variant_a.strandA != structural_variant_a.strandB or \
                structural_variant_b.strandA != structural_variant_b.strandB:
            structural_variant_a = get_simple_sv_from_breakend(breakend_a)
            structural_variants.append(structural_variant_a)
            structural_variant_b = get_simple_sv_from_breakend(breakend_b)
            structural_variants.append(structural_variant_b)
            variant_records.append(variant_record_a)
            variant_records.append(variant_record_b)
        else:
            # put the structural variants all in + orientation, if not already
            structural_variant_a = positive_strand(structural_variant_a)
            structural_variant_b = positive_strand(structural_variant_b)
            # check if the pair of SVs form one dispersed duplication
            if check_dispersed_duplication and \
                    structural_variant_a.chromosomeA == structural_variant_b.chromosomeB and \
                    structural_variant_a.chromosomeB == structural_variant_b.chromosomeA:
                if abs(structural_variant_a.startA - structural_variant_b.startB) <= 10 and \
                        structural_variant_a.startB < structural_variant_b.startA:
                        # dispersed duplication, extract right coordinates
                    chrom_A = structural_variant_a.chromosomeB
                    start_A = structural_variant_a.startB
                    end_A = structural_variant_b.startA
                    chrom_B = structural_variant_a.chromosomeA
                    start_B = structural_variant_a.startA
                    end_B = structural_variant_a.startA
                    structural_variant = SV(
                        chromosomeA=chrom_A,
                        startA=start_A, endA=end_A,
                        sv_type="DUP:DISPERSED",
                        chromosomeB=chrom_B,
                        startB=start_B, endB=end_B)
                    structural_variants.append(structural_variant)
                    # take variant record from a single breakend
                    breakend = breakend_cluster[0]
                    variant_record = breakend.variant_record
                    variant_records.append(variant_record)
                elif abs(structural_variant_a.startB - structural_variant_b.startA) <= 10 and \
                        structural_variant_a.startA > structural_variant_b.startB:
                    # dispersed duplication, extract right coordinates
                    chrom_A = structural_variant_a.chromosomeA
                    start_A = structural_variant_b.startB
                    end_A = structural_variant_a.startA
                    chrom_B = structural_variant_a.chromosomeB
                    start_B = structural_variant_a.startB
                    end_B = structural_variant_a.startB
                    structural_variant = SV(
                        chromosomeA=chrom_A,
                        startA=start_A, endA=end_A,
                        sv_type="DUP:DISPERSED",
                        chromosomeB=chrom_B,
                        startB=start_B, endB=end_B)
                    structural_variants.append(structural_variant)
                    # take variant record from a single breakend
                    breakend = breakend_cluster[0]
                    variant_record = breakend.variant_record
                    variant_records.append(variant_record)
            else:
                # not a dispersed duplication
                structural_variant_a = get_simple_sv_from_breakend(breakend_a)
                structural_variants.append(structural_variant_a)
                structural_variant_b = get_simple_sv_from_breakend(breakend_b)
                structural_variants.append(structural_variant_b)
                variant_records.append(variant_record_a)
                variant_records.append(variant_record_b)
    else:
        # go through SVs one by one, making sure to not process mates
        for breakend in breakend_cluster:
            if breakend.identifier in mate_ids:
                continue
            # store mate id
            mate_ids.append(breakend.mate_id)
            structural_variant = get_simple_sv_from_breakend(breakend)
            structural_variants.append(structural_variant)
            # get variant record
            variant_records.append(breakend.variant_record)
    # add features to each structural variant
    final_structural_variants = []
    for i, sv in enumerate(structural_variants):
        final_sv = add_features(sv, variant_records[i])
        final_structural_variants.append(final_sv)
    return final_structural_variants

def write_structural_variants_tobedpe(structural_variants, bedpe_fn, tool):
    """
    Write structural variant to BEDPE file

    :param structural_variants: List of structural variants
    :param bedpe_fn: Name of output BEDPE file
    :param tool: Tool used to call variants
    :return: 0
    """
    # sort structural variants
    structural_variants.sort(key=lambda x: [x.chromosomeA, x.startA, x.endA])
    with open(bedpe_fn, "w") as output_file:
        # write header to output bed
        header = ["#CHROM_A", "START_A", "END_A", "CHROM_B", "START_B",
                  "END_B",
                  "ID", "QUAL", "STRAND_A", "STRAND_B", "TYPE",
                  "INSERTED_SEQUENCE", "TOOL",
                  "READ_PAIRS_DEFINED",
                  "READ_PAIRS", "SPLIT_READS_DEFINED", "SPLIT_READS",
                  "DELLY_LOW_QUAL_DEFINED",
                  "DELLY_LOW_QUAL", "MANTA_REF_READ_PAIRS_DEFINED",
                  "MANTA_REF_READ_PAIRS", "MANTA_REF_SPLIT_READS_DEFINED",
                  "MANTA_REF_SPLIT_READS"]
        header = '\t'.join(header) + '\n'
        output_file.write(header)
        # write all structural variants to bed file
        for variant in structural_variants:
            if tool == "DELLY":
                variant.delly_low_qual_defined = 1
            # correctly write coordinates in BEDPE format for BND
            if variant.sv_type == "BND":
                if variant.strandA == "-":
                    variant.startA = variant.startA - 1
                    variant.endA = variant.endA - 1
                if variant.strandB == "+":
                    variant.startB = variant.startB - 1
                    variant.endB = variant.endB - 1
            # correctly write coordinates in BEDPE format for INS
            elif variant.sv_type == "INS":
                variant.startA = variant.startA
                variant.endA = variant.endA - 1
            else:
                # convert start position to 0-based coordinate, if not BND or INS
                variant.startA = variant.startA - 1
            output_line_list = [variant.chromosomeA, str(variant.startA),
                                str(variant.endA), variant.chromosomeB,
                                str(variant.startB),
                                str(variant.endB),
                                variant.identifier, str(variant.qual),
                                variant.strandA, variant.strandB,
                                variant.sv_type, variant.inserted_sequence, tool, str(variant.read_pairs_defined),
                                variant.read_pairs, str(variant.split_reads_defined),
                                variant.split_reads, str(variant.delly_low_qual_defined),
                                str(variant.delly_low_qual), str(variant.manta_ref_read_pairs_defined),
                                str(variant.manta_ref_read_pairs), str(variant.manta_ref_split_reads_defined),
                                str(variant.manta_ref_split_reads)]
            output_line = '\t'.join(output_line_list) + '\n'
            output_file.write(output_line)
    return 0

def process_vcf(input_vcf_fn, output_bedpe_fn, tool, check_dispersed_duplication=True):
    """Process vcfs, extracting simple structural variants

    :param input_vcf_fn: Filename of input VCF file
    :param output_vcf_fn: Filename of output BEDPE file
    :param tool: Tool used to call variants
    :param check_dispersed_duplication: Will check for dispersed duplications
    if True (boolean)
    :return: 0 (integer)
    """
    # collect breakends and insertions from vcf file
    with pysam.VariantFile(input_vcf_fn) as vcf:
        breakends = parse_breakends_and_insertions_from_vcf(vcf)
        # loop through each line of the vcf
    # cluster breakends
    breakend_clusters = cluster_breakends(breakends)
    # extract structural variants from clusters
    final_structural_variants = []
    for breakend_cluster in breakend_clusters.values():
        structural_variants = structural_variant_from_breakend_cluster(breakend_cluster, check_dispersed_duplication)
        if structural_variants:
            for variant in structural_variants:
                final_structural_variants.append(variant)
    # write structural variants to bedpe file
    write_structural_variants_tobedpe(final_structural_variants, output_bedpe_fn, tool)
    return 0

def main():
    args = parse_cl_args(sys.argv[1:])
    # convert vcf file
    process_vcf(args.input_vcf, args.output_bedpe, args.tool,
                args.dispersed_duplications)

if __name__ == "__main__":
    main()
