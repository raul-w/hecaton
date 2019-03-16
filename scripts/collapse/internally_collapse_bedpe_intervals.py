#!/usr/bin/env python3

"""
Collapse all overlapping variants within a BEDPE file
"""

import argparse
import pandas as pd
import sys
from intervaltree import Interval, IntervalTree
from utils.utils import get_overlap
from itertools import count
from classes.sv import SV
from statistics import mean
from math import ceil


def parse_cl_args(in_args):
    """
    Parse command line arguments

    :param in_args: All command line arguments
    :return: None
    """
    description = "Collapse all overlapping variants within a BEDPE file"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-b", "--bedpe_fn", type=str,
                        help="Path to BEDPE file B")
    parser.add_argument("-f", "--fai_fn", type=str,
                        help="Path to fasta index file of used genome")
    parser.add_argument("-o", "--output_file", type=str,
                        help="Output file containing intersected intervals")
    parser.add_argument("-e", "--extra_features", type=bool, default=False,
                        help="Indicates whether BEDPE files contain "
                             "'extra features'")
    parser.add_argument("-r", "--reciprocal", type=float, default=0.5,
                        help="Minimum fraction of reciprocal overlap needed "
                             "to collapse calls")
    args = parser.parse_args(in_args)
    return args

def check_overlap(cnv, svs, chrom_length_dict, id_generator, extra_features=False, reciprocal=0.5):
    """Check if CNV is present in SV dictionary, update entry if it is, add
    new entry if it is not

    :param cnv: SV instance
    :param svs: Dictionary containing chromosome names as keys and
    interval trees of SVs as values
    :param chrom_length_dict: Dictionary with chromosome names as keys and
    their lengths as values
    :param id_generator: Count object from itertools
    :param extra_features: Boolean indicating whether CNVs contain 'extra
    features'
    :param reciprocal: Minimum fraction of reciprocal overlap needed to collapse calls (float)
    :return: Updated version of svs, updated count object
    """
    # skip BNDs
    if cnv.TYPE == "BND":
        return svs, id_generator
    chromA = str(cnv[0])
    startA = cnv.START_A + 1
    endA = cnv.END_A
    chromB = cnv.CHROM_B
    startB = cnv.START_B
    endB = cnv.END_B
    qual = cnv.QUAL
    strand_a = cnv.STRAND_A
    strand_b = cnv.STRAND_B
    sv_type = cnv.TYPE
    inserted_sequence = cnv.INSERTED_SEQUENCE
    tool = cnv.TOOL
    read_pairs_defined = cnv.READ_PAIRS_DEFINED
    read_pairs = str(cnv.READ_PAIRS)
    split_reads_defined = cnv.SPLIT_READS_DEFINED
    split_reads = str(cnv.SPLIT_READS)
    if extra_features:
        delly_low_qual_defined = cnv.DELLY_LOW_QUAL_DEFINED
        delly_low_qual = cnv.DELLY_LOW_QUAL
        manta_ref_read_pairs_defined = cnv.MANTA_REF_READ_PAIRS_DEFINED
        manta_ref_read_pairs = cnv.MANTA_REF_READ_PAIRS
        manta_ref_split_reads_defined = cnv.MANTA_REF_SPLIT_READS_DEFINED
        manta_ref_split_reads = cnv.MANTA_REF_SPLIT_READS
    # check if chromosome is present in interval tree
    if chromA not in svs:
        svs[chromA] = IntervalTree()
    # check if cnv has overlap with any previously detected ones
    cnv_interval = [startA, endA]
    cnv_interval_length = endA - startA + 1
    interval_start = cnv_interval[0] - 1
    interval_end = cnv_interval[1]
    # check overlap with point in case of insertion
    if interval_start == interval_end:
        overlapping_cnvs = svs[chromA].at(interval_start)
    else:
        # get calls that overlap with the midpoint or a slightly wider range
        # aim is to reduce the search range as much as possible to speed up queries
        # get midpoint
        if cnv_interval_length % 2 == 0:
            midpoint_start = mean([interval_start, interval_end])
            midpoint_end = mean([interval_start, interval_end])
        else:
            midpoint_start = mean([interval_start, interval_end]) - 0.5
            midpoint_end = mean([interval_start, interval_end]) + 0.5
        reciprocal_slop = ceil(reciprocal * cnv_interval_length)
        reciprocal_interval_start = min(interval_start + reciprocal_slop, interval_end - reciprocal_slop)
        reciprocal_interval_end = max(interval_start + reciprocal_slop, interval_end - reciprocal_slop)
        # search smallest range
        if midpoint_start >= reciprocal_interval_start:
            if midpoint_start == midpoint_end:
                overlapping_cnvs = svs[chromA].at(midpoint_start)
            else:
                overlapping_cnvs = svs[chromA].overlap(midpoint_start, midpoint_end)
        else:
            overlapping_cnvs = svs[chromA].overlap(reciprocal_interval_start, reciprocal_interval_end)
    # check if there are any cnvs of the same type, saving only the best hit
    best_hit = None
    for old_cnv_interval_instance in overlapping_cnvs:
        overlap_found = False
        old_cnv = old_cnv_interval_instance.data
        # compare types
        duplication_types = ["DUP:TANDEM", "DUP:DISPERSED", "DUP"]
        if sv_type in duplication_types and old_cnv.sv_type in duplication_types:
            # duplications are a special case, only a tandem duplication
            # and dispersed duplication will not match
            if sv_type == "DUP:TANDEM" and old_cnv.sv_type == "DUP:DISPERSED":
                continue
            elif sv_type == "DUP:DISPERSED" and old_cnv.sv_type == "DUP:TANDEM":
                continue
            # skip dispersed duplications whose insertion sites differ by more than 10 bp
            if sv_type == "DUP:DISPERSED" and old_cnv.sv_type == "DUP:DISPERSED":
                if abs(startB - old_cnv.startB) >= 10:
                    continue
        else:
            # not a duplication, just compare types
            if sv_type != old_cnv.sv_type:
                continue
        if sv_type == "INS":
            # insertions must be within 10 bp if there is overlap
            overlap_found = True
        else:
            # check if there is enough reciprocal overlap and breakpoints are within 1 kbp
            old_cnv_interval = [old_cnv.startA, old_cnv.endA]
            old_cnv_interval_length = old_cnv_interval[1] - old_cnv_interval[
                0] + 1
            overlap = get_overlap(old_cnv_interval, cnv_interval)
            if overlap > reciprocal * cnv_interval_length and overlap > reciprocal * old_cnv_interval_length \
                    and abs(cnv_interval[0] - old_cnv_interval[0]) <= 1000 \
                    and abs(cnv_interval[1] - old_cnv_interval[1]) <= 1000:
                overlap_found = True
        if overlap_found:
            if best_hit is None:
                best_hit = old_cnv_interval_instance
            else:
                # compare support of this hit with that of the best hit and replace best hit if necessary
                best_hit_data = best_hit.data
                if old_cnv.split_reads > best_hit_data.split_reads:
                    best_hit = old_cnv_interval_instance
                elif old_cnv.split_reads < best_hit_data.split_reads:
                    continue
                elif old_cnv.read_pairs > best_hit_data.read_pairs:
                    best_hit = old_cnv_interval_instance
                elif old_cnv.read_pairs < best_hit_data.read_pairs:
                    continue
                else:
                    # support is tied, just keep the previously found best hit
                    continue
    if best_hit:
        # update entry in svs dictionary
        best_hit_cnv = best_hit.data
        new_cnv = best_hit_cnv
        # take coordinates of best supported interval
        best_hit_most_supported = False
        if split_reads > best_hit_cnv.split_reads:
            pass
        elif split_reads < best_hit_cnv.split_reads:
            best_hit_most_supported = True
        elif read_pairs > best_hit_cnv.read_pairs:
            pass
        elif read_pairs < best_hit_cnv.read_pairs:
            best_hit_most_supported = True
        else:
            # a tie, take the coordinates of the best hit
            best_hit_most_supported = True
        if not best_hit_most_supported:
            # replace coordinates of new cnv
            new_cnv.startA = startA
            # ensure end does not exceed chrom length
            new_cnv.endA = min(endA, chrom_length_dict[str(chromA)])
            # add read pairs and split reads
            new_cnv.read_pairs_defined = read_pairs_defined
            new_cnv.read_pairs = read_pairs
            new_cnv.split_reads_defined = split_reads_defined
            new_cnv.split_reads = split_reads
            if extra_features:
                new_cnv.delly_low_qual_defined = delly_low_qual_defined
                new_cnv.delly_low_qual = delly_low_qual
                new_cnv.manta_ref_read_pairs_defined = manta_ref_read_pairs_defined
                new_cnv.manta_ref_read_pairs = manta_ref_read_pairs
                new_cnv.manta_ref_split_reads_defined = manta_ref_split_reads_defined
                new_cnv.manta_ref_split_reads = manta_ref_split_reads
        new_tools = ';'.join([best_hit_cnv.tool, tool])
        new_cnv.tool = new_tools
        # make DUP more specific, if possible
        if new_cnv.sv_type == "DUP":
            new_cnv.sv_type = sv_type
        # remove old cnv interval
        svs[chromA].remove(best_hit)
        # add new interval
        if new_cnv.sv_type == "INS":
            interval_start = max(0, startA - 11)
            interval_end = min(endA + 10, chrom_length_dict[str(chromA)])
            svs[chromA].addi(interval_start, interval_end, new_cnv)
        else:
            svs[chromA].addi(new_cnv.startA - 1, new_cnv.endA, new_cnv)
    else:
        # add new cnv
        identifier = next(id_generator)
        # ensure that end does not exceed chromosome
        endA = min(endA, chrom_length_dict[str(chromA)])
        new_sv = SV(chromosomeA=chromA, startA=startA, endA=endA,
                    chromosomeB=chromB, startB=startB, endB=endB,
                    identifier=identifier, qual=qual, strandA=strand_a,
                    strandB=strand_b, sv_type=sv_type,
                    inserted_sequence=inserted_sequence,
                    tool=tool,
                    read_pairs_defined=read_pairs_defined,
                    read_pairs=read_pairs,
                    split_reads_defined=split_reads_defined,
                    split_reads=split_reads)
        if extra_features:
            new_sv.delly_low_qual_defined = delly_low_qual_defined
            new_sv.delly_low_qual = delly_low_qual
            new_sv.manta_ref_read_pairs_defined = manta_ref_read_pairs_defined
            new_sv.manta_ref_read_pairs = manta_ref_read_pairs
            new_sv.manta_ref_split_reads_defined = manta_ref_split_reads_defined
            new_sv.manta_ref_split_reads = manta_ref_split_reads
        # add new interval
        if new_sv.sv_type == "INS":
            # add 10 bp upstream and downstream for overlap queries
            interval_start = max(0, startA - 11)
            interval_end = min(endA + 10, chrom_length_dict[str(chromA)])
            svs[chromA].addi(interval_start, interval_end, new_sv)
        else:
            svs[chromA].addi(new_sv.startA - 1, new_sv.endA, new_sv)
    return svs, id_generator

def collapse_bedpe(bedpe_fn, output_bedpe_fn, chrom_length_dict, extra_features, reciprocal=0.5):
    """
    Collapse all overlapping variants within a BEDPE file

    :param bedpe_fn: Filename of input BEDPE file
    :param output_bedpe_fn: Filename of output BEDPE file
    :param chrom_length_dict: Dictionary with chromosomes as keys and
    length of chromosomes as values
    :param extra_features: Boolean indicating whether CNVs contain 'extra
    features'
    :param reciprocal: Minimum fraction of reciprocal overlap needed to collapse calls (float)
    :return: 0 (integer)
    """
    # get svs out of bedpe file A, giving them new identifiers
    id_generator = count(0)
    svs = {}
    bedpe_a = pd.read_csv(bedpe_fn, sep="\t")
    for cnv in bedpe_a.itertuples(index=False, name='Pandas'):
        svs, id_generator = check_overlap(cnv, svs, chrom_length_dict, id_generator, extra_features, reciprocal)
    # write all svs to bedpe
    structural_variants = []
    for chrom in svs:
        for interval in svs[chrom]:
            structural_variants.append(interval.data)
    # sort structural variants, if not empty
    if structural_variants:
        structural_variants.sort(
            key=lambda x: [x.chromosomeA, x.startA, x.endA])
    with open(output_bedpe_fn, "w") as output_file:
        # write header to output bed
        if extra_features:
            header = ["#CHROM_A", "START_A", "END_A", "CHROM_B", "START_B",
                      "END_B",
                      "ID", "QUAL", "STRAND_A", "STRAND_B", "TYPE",
                      "INSERTED_SEQUENCE", "TOOL",
                      "READ_PAIRS_DEFINED",
                      "READ_PAIRS", "SPLIT_READS_DEFINED", "SPLIT_READS",
                      "DELLY_LOW_QUAL_DEFINED",
                      "DELLY_LOW_QUAL", "MANTA_REF_READ_PAIRS_DEFINED",
                      "MANTA_REF_READ_PAIRS", "MANTA_REF_SPLIT_READS_DEFINED",
                      "MANTA_REF_SPLIT_READS", "SIZE"]
        else:
            header = ["#CHROM_A", "START_A", "END_A", "CHROM_B", "START_B",
                  "END_B",
                  "ID", "QUAL", "STRAND_A", "STRAND_B", "TYPE",
                  "INSERTED_SEQUENCE", "TOOL",
                  "READ_PAIRS_DEFINED", "READ_PAIRS", "SPLIT_READS_DEFINED",
                  "SPLIT_READS", "SIZE"]
        header = '\t'.join(header) + '\n'
        output_file.write(header)
        # write all structural variants to bed file
        for variant in structural_variants:
            # create correct size
            if variant.sv_type == "INS":
                if variant.inserted_sequence == "N":
                    size = 5000
                else:
                    size = len(variant.inserted_sequence)
            else:
                size = variant.endA - variant.startA + 1
            # convert start position to 0-based coordinate
            if extra_features:
                output_line_list = [variant.chromosomeA,
                                    str(variant.startA - 1),
                                    str(variant.endA), variant.chromosomeB,
                                    str(variant.startB), str(variant.endB),
                                    str(variant.identifier), str(variant.qual),
                                    variant.strandA, variant.strandB,
                                    variant.sv_type, variant.inserted_sequence,
                                    variant.tool,
                                    str(variant.read_pairs_defined),
                                    str(variant.read_pairs),
                                    str(variant.split_reads_defined),
                                    str(variant.split_reads),
                                    str(variant.delly_low_qual_defined),
                                    str(variant.delly_low_qual),
                                    str(variant.manta_ref_read_pairs_defined),
                                    str(variant.manta_ref_read_pairs),
                                    str(variant.manta_ref_split_reads_defined),
                                    str(variant.manta_ref_split_reads),
                                    str(size)]
            else:
                output_line_list = [variant.chromosomeA, str(variant.startA - 1),
                                    str(variant.endA), variant.chromosomeB,
                                    str(variant.startB), str(variant.endB),
                                    str(variant.identifier), str(variant.qual),
                                    variant.strandA, variant.strandB,
                                    variant.sv_type, variant.inserted_sequence,
                                    variant.tool,
                                    str(variant.read_pairs_defined),
                                    str(variant.read_pairs),
                                    str(variant.split_reads_defined),
                                    str(variant.split_reads), str(size)]
            output_line = '\t'.join(output_line_list) + '\n'
            output_file.write(output_line)
    return 0

def main():
    args = parse_cl_args(sys.argv[1:])
    # check if reciprocal is between 0 and 1
    if args.reciprocal < 0 or args.reciprocal > 1:
        raise ValueError("Reciprocal overlap should be fraction between 0 and 1")
    # get length of all chromosomes
    chrom_length_dict = {}
    with open(args.fai_fn) as fai_file:
        for line in fai_file:
            line_elems = line.strip().split()
            chrom = line_elems[0]
            chrom_length = int(line_elems[1])
            chrom_length_dict[chrom] = chrom_length
    # collapse bedpe
    collapse_bedpe(args.bedpe_fn, args.output_file,
                    chrom_length_dict, args.extra_features, args.reciprocal)

if __name__ == "__main__":
    main()