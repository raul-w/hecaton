#!/usr/bin/env python3

"""
Intersect two BEDPE files that contain SVs
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
    description = "Parse VCF files, extracting simple CNVs in BEDPE format"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-a", "--bedpe_a", type=str,
                        help="Path to BEDPE file A")
    parser.add_argument("-b", "--bedpe_b", type=str,
                        help="Path to BEDPE file B")
    parser.add_argument("-f", "--fai_fn", type=str,
                        help="Path to fasta index file of used genome")
    parser.add_argument("-o", "--output_file", type=str,
                        help="Output file containing intersected intervals")
    parser.add_argument("-r", "--reciprocal", type=float, default=0.5,
                        help="Minimum fraction of reciprocal overlap needed "
                             "to collapse calls")
    args = parser.parse_args(in_args)
    return args

def intersect_bedpe(bedpe_a_fn, bedpe_b_fn, output_bedpe_fn, chrom_length_dict, reciprocal=0.5):
    """
    Intersect two BEDPE files using interval tree

    :param bedpe_a_fn: Filename of input BEDPE file A
    :param bedpe_b_fn: Filename of input BEDPE file B
    :param output_bedpe_fn: Filename of output BEDPE file
    :param chrom_length_dict: Dictionary with chromosomes as keys and
    length of chromosomes as values
    :return: 0 (integer)
    """
    # get svs out of bedpe file A, giving them new identifiers
    id_generator = count(0)
    svs = {}
    bedpe_a = pd.read_csv(bedpe_a_fn, sep="\t")
    for cnv in bedpe_a.itertuples(index=False, name='Pandas'):
        # skip BNDs
        if cnv.TYPE == "BND":
            continue
        chromA = str(cnv[0])
        startA = cnv.START_A + 1
        endA = cnv.END_A
        chromB = cnv.CHROM_B
        startB = cnv.START_B
        endB = cnv.END_B
        identifier = next(id_generator)
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
        delly_low_qual_defined = cnv.DELLY_LOW_QUAL_DEFINED
        delly_low_qual = cnv.DELLY_LOW_QUAL
        manta_ref_read_pairs_defined = cnv.MANTA_REF_READ_PAIRS_DEFINED
        manta_ref_read_pairs = cnv.MANTA_REF_READ_PAIRS
        manta_ref_split_reads_defined = cnv.MANTA_REF_SPLIT_READS_DEFINED
        manta_ref_split_reads = cnv.MANTA_REF_SPLIT_READS
        # check if chromosome is present in interval tree
        if chromA not in svs:
            svs[chromA] = IntervalTree()
        new_sv = SV(chromosomeA=chromA, startA=startA, endA=endA,
                    chromosomeB=chromB, startB=startB, endB=endB,
                    identifier=identifier, qual=qual, strandA=strand_a,
                    strandB=strand_b, sv_type=sv_type, inserted_sequence=inserted_sequence,
                    tool=tool, read_pairs_defined=read_pairs_defined,
                    read_pairs=read_pairs, split_reads_defined=split_reads_defined, split_reads=split_reads,
                    delly_low_qual_defined=delly_low_qual_defined, delly_low_qual=delly_low_qual,
                    manta_ref_read_pairs_defined=manta_ref_read_pairs_defined, manta_ref_read_pairs=manta_ref_read_pairs,
                    manta_ref_split_reads_defined=manta_ref_split_reads_defined, manta_ref_split_reads=manta_ref_split_reads)
        if new_sv.sv_type == "INS":
            # add 10 bp upstream and downstream for overlap queries
            interval_start = max(0, startA - 11)
            interval_end = min(endA + 10, chrom_length_dict[str(chromA)])
            svs[chromA].addi(interval_start, interval_end, new_sv)
        else:
            svs[chromA].addi(new_sv.startA - 1, new_sv.endA, new_sv)
    # search for more than 50% reciprocal overlap between bedpe a and b
    bedpe_b = pd.read_csv(bedpe_b_fn, sep="\t")
    for cnv in bedpe_b.itertuples(index=False, name='Pandas'):
        # skip BNDs
        if cnv.TYPE == "BND":
            continue
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
        read_pairs_defined = int(cnv.READ_PAIRS_DEFINED)
        read_pairs = str(cnv.READ_PAIRS)
        split_reads_defined = int(cnv.SPLIT_READS_DEFINED)
        split_reads = str(cnv.SPLIT_READS)
        delly_low_qual_defined = int(cnv.DELLY_LOW_QUAL_DEFINED)
        delly_low_qual = int(cnv.DELLY_LOW_QUAL)
        manta_ref_read_pairs_defined = int(cnv.MANTA_REF_READ_PAIRS_DEFINED)
        manta_ref_read_pairs = int(cnv.MANTA_REF_READ_PAIRS)
        manta_ref_split_reads_defined = int(cnv.MANTA_REF_SPLIT_READS_DEFINED)
        manta_ref_split_reads = int(cnv.MANTA_REF_SPLIT_READS)
        # check if chromosome is present in interval tree
        if chromA not in svs:
            svs[chromA] = IntervalTree()
        # check if there are any cnvs with enough reciprocal overlap
        cnv_interval = [startA, endA]
        cnv_interval_length = endA - startA + 1
        interval_start = cnv_interval[0] - 1
        interval_end = cnv_interval[1]
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
            reciprocal_interval_start = min(interval_start + reciprocal_slop,
                                            interval_end - reciprocal_slop)
            reciprocal_interval_end = max(interval_start + reciprocal_slop,
                                          interval_end - reciprocal_slop)
            # search smallest range
            if midpoint_start >= reciprocal_interval_start:
                if midpoint_start == midpoint_end:
                    overlapping_cnvs = svs[chromA].at(midpoint_start)
                else:
                    overlapping_cnvs = svs[chromA].overlap(midpoint_start,
                                                           midpoint_end)
            else:
                overlapping_cnvs = svs[chromA].overlap(
                    reciprocal_interval_start, reciprocal_interval_end)
        # check if there are any cnvs of the same type with enough reciprocal overlap, saving only the best hit
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
                # skip dispersed duplications whose insertion sites differ by more than 20 bp
                if sv_type == "DUP:DISPERSED" and old_cnv.sv_type == "DUP:DISPERSED":
                    if abs(startB - old_cnv.startB) >= 10:
                        continue
            else:
                # not a duplication, just compare types
                if sv_type != old_cnv.sv_type:
                    continue
            old_cnv_interval = [old_cnv.startA, old_cnv.endA]
            # insertions need specific treatment
            if sv_type == "INS":
                # insertions must be within 10 bp if there is overlap
                overlap_found = True
            else:
                # check if there is more than 50% reciprocal overlap and breakpoints are within 1 kbp
                old_cnv_interval_length = old_cnv_interval[1] - old_cnv_interval[0] + 1
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
            # update coordinates, taking the union of the intervals
            if sv_type == "INS":
                new_cnv.startA = min(best_hit_cnv.startA, startA)
                new_cnv.endA = max(best_hit_cnv.endA, endA)
            else:
                new_cnv.startA = min(best_hit_cnv.startA, startA)
                new_cnv.endA = max(best_hit_cnv.endA, endA)
                # ensure end does not exceed chrom length
                new_cnv.endA = min(best_hit_cnv.endA,
                                   chrom_length_dict[str(chromA)])
            new_tools = ';'.join([best_hit_cnv.tool, tool])
            new_cnv.tool = new_tools
            # make DUP more specific, if possible
            if new_cnv.sv_type == "DUP":
                new_cnv.sv_type = sv_type
            # join split reads and read pairs
            new_cnv.read_pairs_defined = max(read_pairs_defined, best_hit_cnv.read_pairs_defined)
            new_cnv.read_pairs = ";".join([best_hit_cnv.read_pairs, read_pairs])
            new_cnv.split_reads_defined = max(split_reads_defined, best_hit_cnv.split_reads_defined)
            new_cnv.split_reads = ";".join([best_hit_cnv.split_reads, split_reads])
            new_cnv.delly_low_qual_defined = max(delly_low_qual_defined,
                                                 best_hit_cnv.delly_low_qual_defined)
            new_cnv.delly_low_qual = max(delly_low_qual, best_hit_cnv.delly_low_qual)
            new_cnv.manta_ref_read_pairs_defined = max(manta_ref_read_pairs_defined,
                                                       best_hit_cnv.manta_ref_read_pairs_defined)
            new_cnv.manta_ref_read_pairs = max(manta_ref_read_pairs, best_hit_cnv.manta_ref_read_pairs)
            new_cnv.manta_ref_split_reads_defined = max(manta_ref_split_reads_defined,
                                                        best_hit_cnv.manta_ref_split_reads_defined)
            new_cnv.manta_ref_split_reads = max(manta_ref_split_reads,
                                                best_hit_cnv.manta_ref_split_reads)
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
            identifier = next(id_generator)
            # ensure that end does not exceed chromosome
            endA = min(endA, chrom_length_dict[str(chromA)])
            new_sv = SV(chromosomeA=chromA, startA=startA, endA=endA,
                    chromosomeB=chromB, startB=startB, endB=endB,
                    identifier=identifier, qual=qual, strandA=strand_a,
                    strandB=strand_b, sv_type=sv_type, inserted_sequence=inserted_sequence,
                    tool=tool, read_pairs_defined=read_pairs_defined,
                    read_pairs=read_pairs, split_reads_defined=split_reads_defined, split_reads=split_reads,
                    delly_low_qual_defined=delly_low_qual_defined, delly_low_qual=delly_low_qual,
                    manta_ref_read_pairs_defined=manta_ref_read_pairs_defined, manta_ref_read_pairs=manta_ref_read_pairs,
                    manta_ref_split_reads_defined=manta_ref_split_reads_defined, manta_ref_split_reads=manta_ref_split_reads)
            # add new interval
            if new_sv.sv_type == "INS":
                # add 10 bp upstream and downstream for overlap queries
                interval_start = max(0, startA - 11)
                interval_end = min(endA + 10, chrom_length_dict[str(chromA)])
                svs[chromA].addi(interval_start, interval_end, new_sv)
            else:
                svs[chromA].addi(new_sv.startA - 1, new_sv.endA, new_sv)
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
        header = ["#CHROM_A", "START_A", "END_A", "CHROM_B", "START_B",
                  "END_B",
                  "ID", "QUAL", "STRAND_A", "STRAND_B", "TYPE",
                  "INSERTED_SEQUENCE", "TOOL", "READ_PAIRS_DEFINED",
                  "READ_PAIRS", "SPLIT_READS_DEFINED", "SPLIT_READS",
                  "DELLY_LOW_QUAL_DEFINED",
                  "DELLY_LOW_QUAL", "MANTA_REF_READ_PAIRS_DEFINED",
                  "MANTA_REF_READ_PAIRS", "MANTA_REF_SPLIT_READS_DEFINED",
                  "MANTA_REF_SPLIT_READS"]
        header = '\t'.join(header) + '\n'
        output_file.write(header)
        # write all structural variants to bed file
        for variant in structural_variants:
            # convert start position to 0-based coordinate
            output_line_list = [variant.chromosomeA, str(variant.startA - 1),
                                str(variant.endA), variant.chromosomeB,
                                str(variant.startB), str(variant.endB),
                                str(variant.identifier), str(variant.qual),
                                variant.strandA, variant.strandB,
                                variant.sv_type, variant.inserted_sequence,
                                variant.tool, str(variant.read_pairs_defined),
                                str(variant.read_pairs), str(variant.split_reads_defined),
                                str(variant.split_reads), str(variant.delly_low_qual_defined),
                                str(variant.delly_low_qual), str(variant.manta_ref_read_pairs_defined),
                                str(variant.manta_ref_read_pairs), str(variant.manta_ref_split_reads_defined),
                                str(variant.manta_ref_split_reads)]
            output_line = '\t'.join(output_line_list) + '\n'
            output_file.write(output_line)
    return 0

def main():
    args = parse_cl_args(sys.argv[1:])
    # get length of all chromosomes
    chrom_length_dict = {}
    with open(args.fai_fn) as fai_file:
        for line in fai_file:
            line_elems = line.strip().split()
            chrom = line_elems[0]
            chrom_length = int(line_elems[1])
            chrom_length_dict[chrom] = chrom_length
    # convert vcf file
    intersect_bedpe(args.bedpe_a, args.bedpe_b, args.output_file,
                    chrom_length_dict)

if __name__ == "__main__":
    main()