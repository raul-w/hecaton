#!/usr/bin/env python3

"""
Merge VCF files of single samples containing CNVs called by Hecaton, producing
one merged VCF file
"""

import argparse
import datetime
import operator
import pysam
import sys
from collections import defaultdict
from intervaltree import Interval, IntervalTree
from itertools import count
from statistics import mean
from math import ceil
from utils.utils import get_overlap


def parse_cl_args(in_args):
    """
    Parse command line arguments

    :param in_args: All command line arguments
    :return: None
    """
    description = "Merge VCF files of single samples containing CNVs called by " \
                  "Hecaton, producing VCF file containing the merged calls"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-i", "--input_file", type=str,
                        help="File containing path to single sample VCF file "
                             "on each line")
    parser.add_argument("-f", "--fai_fn", type=str,
                        help="Path to fasta index file of used genome")
    parser.add_argument("-o", "--output_vcf", type=str,
                        help="Name of VCF output file")
    parser.add_argument("-r", "--reciprocal", type=float, default=0.5,
                        help="Minimum fraction of reciprocal overlap needed "
                             "to collapse calls")
    args = parser.parse_args(in_args)
    return args

def obtain_sites_and_genotypes(input_fns):
    """
    Obtain all CNV sites from a list of VCF files

    :param input_fns: List of paths to VCF files
    :return: Set containing all sites. Sites are tuples
    (chrom, pos, alt, end, sv_type, inschrom, inspos). Dict containing sites
    as keys and genotype lines of samples as values.
    """
    # obtain all sites
    sites_set = set()
    # obtain all genotypes of samples
    samples_dict = {}
    for input_vcf_fn in input_fns:
        with pysam.VariantFile(input_vcf_fn) as vcf:
            # add new sample to sample dict
            sample = vcf.header.samples[0]
            samples_dict[sample] = defaultdict(int)
            for record in vcf.fetch():
                # get site information
                chrom = record.chrom
                pos = record.pos
                alt = record.alts[0]
                end = record.stop
                sv_type = record.info["SVTYPE"]
                # check for dispersed duplications
                if sv_type == "DUP:DISPERSED":
                    inschrom = record.info["INSCHROM"][0]
                    inspos = record.info["INSPOS"][0]
                else:
                    inschrom = ""
                    inspos = -1
                site_tuple = (chrom, pos, alt, end, sv_type, inschrom, inspos)
                sites_set.add(site_tuple)
                # get genotype information
                sup = str(record.samples[0]["SUP"])
                rp = str(record.samples[0]["RP"])
                sr = str(record.samples[0]["SR"])
                tool = ""
                if "TOOL" in record.format:
                    tool = record.samples[0]["TOOL"]
                    tool = ",".join(list(tool))
                rq = str(round(record.samples[0]["RQ"], 2))
                dhfc = record.samples[0]["DHFC"]
                dhbfc = record.samples[0]["DHBFC"]
                dhffc = record.samples[0]["DHFFC"]
                # set genotypes based on format
                gt = [".", "."]
                if sv_type == "INS":
                    gt[1] = "1"
                elif sv_type == "DEL":
                    if dhfc >= 0 and dhfc < 0.25:
                        gt[0] = "1"
                        gt[1] = "1"
                    elif dhfc >= 0.25 and dhfc < 0.75:
                        gt[0] = "0"
                        gt[1] = "1"
                    else:
                        gt[0] = "0"
                        gt[1] = "0"
                elif sv_type == "DUP:TANDEM" or sv_type == "DUP:DISPERSED":
                    if dhfc >= 0 and dhfc < 1.25:
                        gt[0] = "0"
                        gt[1] = "0"
                    elif dhfc >= 1.25 and dhfc < 1.75:
                        gt[0] = "0"
                        gt[1] = "1"
                    else:
                        gt[0] = "1"
                        gt[1] = "1"
                else:
                    raise ValueError("Cannot genotype unknown SV type: {}".format(sv_type))
                gt = "/".join(gt)
                # add genotype line to samples dict
                if tool:
                    genotype_line = ":".join([gt, sup, rp, sr, tool, rq, str(dhfc), str(dhbfc), str(dhffc)])
                else:
                    genotype_line = ":".join([gt, sup, rp, sr, rq, str(dhfc), str(dhbfc), str(dhffc)])
                # add genotype information to samples dict
                samples_dict[sample][site_tuple] = genotype_line
    return sites_set, samples_dict

def get_interval_tree_from_sites_set(sites_set, chrom_length_dict, reciprocal=0.5):
    """
    Get interval tree in which sites considered to be the same CNV are
    collapsed

    :param sites_set: Set containing all sites. Sites are tuples
    (chrom, pos, alt, end, sv_type, inschrom, inspos).
    :param chrom_length_dict: Dictionary with chromosomes as keys and
    length of chromosomes as values
    :param reciprocal: Minimum reciprocal overlap at which CNVs are considered
    to be the same variant
    :return: Dictionary containing chromosomes as keys and
    IntervalTree instances containing collapsed sites as values.
    IntervalTree data instances contain dicts with collapsed sites as the key
    and a set of original sites as values.
    """
    sites_interval_trees = {}
    # search for more than 50% reciprocal overlap between bedpe a and b
    for site in sites_set:
        chrom = site[0]
        pos = site[1]
        alt = site[2]
        end = site[3]
        sv_type = site[4]
        inschrom = site[5]
        inspos = site[6]
        # check if chromosome is present in interval tree
        if chrom not in sites_interval_trees:
            sites_interval_trees[chrom] = IntervalTree()
        # check if there are any cnvs with enough reciprocal overlap
        cnv_interval = [pos, end]
        cnv_interval_length = end - pos + 1
        interval_start = cnv_interval[0] - 1
        interval_end = cnv_interval[1]
        if interval_start == interval_end:
            overlapping_cnvs = sites_interval_trees[chrom].at(interval_start)
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
                    overlapping_cnvs = sites_interval_trees[chrom].at(midpoint_start)
                else:
                    overlapping_cnvs = sites_interval_trees[chrom].overlap(midpoint_start,
                                                           midpoint_end)
            else:
                overlapping_cnvs = sites_interval_trees[chrom].overlap(
                    reciprocal_interval_start, reciprocal_interval_end)
        # check if there are any cnvs of the same type with enough reciprocal overlap, saving only the best hit
        best_hit = None
        for old_cnv_interval_instance in overlapping_cnvs:
            overlap_found = False
            old_cnv = next(iter(old_cnv_interval_instance.data))
            # compare types
            duplication_types = ["DUP:TANDEM", "DUP:DISPERSED", "DUP"]
            if sv_type in duplication_types and old_cnv[4] in duplication_types:
                # duplications are a special case, only a tandem duplication
                # and dispersed duplication will not match
                if sv_type == "DUP:TANDEM" and old_cnv[4] == "DUP:DISPERSED":
                    continue
                elif sv_type == "DUP:DISPERSED" and old_cnv[4] == "DUP:TANDEM":
                    continue
                # skip dispersed duplications whose insertion sites differ by more than 10 bp
                if sv_type == "DUP:DISPERSED" and old_cnv[4] == "DUP:DISPERSED":
                    if inschrom != old_cnv[5] or abs(inspos - old_cnv[6]) >= 10:
                        continue
            else:
                # not a duplication, just compare types
                if sv_type != old_cnv[4]:
                    continue
            old_cnv_interval = [old_cnv[1], old_cnv[3]]
            # insertions need specific treatment
            if sv_type == "INS":
                # check if overlap has the same insertion
                if old_cnv[2] == alt:
                    overlap_found = True
            else:
                # check if there is enough reciprocal overlap and breakpoints are within 1 kbp
                old_cnv_interval_length = old_cnv_interval[1] - \
                                          old_cnv_interval[0] + 1
                overlap = get_overlap(old_cnv_interval, cnv_interval)
                if overlap > reciprocal * cnv_interval_length and overlap > reciprocal * old_cnv_interval_length \
                        and abs(cnv_interval[0] - old_cnv_interval[0]) <= 1000 \
                        and abs(cnv_interval[1] - old_cnv_interval[1]) <= 1000:
                    overlap_found = True
            if overlap_found:
                if best_hit is None:
                    best_hit = old_cnv_interval_instance
                else:
                    # compare overlap of this hit with that of the best hit and replace best hit if necessary
                    best_hit_data = next(iter(best_hit.data))
                    best_hit_interval = [best_hit_data[1], best_hit_data[3]]
                    best_hit_overlap = get_overlap(best_hit_interval, cnv_interval)
                    if overlap > best_hit_overlap:
                        best_hit = old_cnv_interval_instance
                    else:
                        # support is less or tied, just keep the previously found best hit
                        continue
        if best_hit:
            # update entry in svs dictionary
            best_hit_cnv = next(iter(best_hit.data))
            new_cnv_original_sites = best_hit.data[best_hit_cnv]
            new_cnv_site = list(best_hit_cnv)
            # update coordinates, taking the union of the intervals
            if sv_type == "INS":
                new_cnv_site[1] = min(best_hit_cnv[1], pos)
                new_cnv_site[3] = max(best_hit_cnv[3], end)
            else:
                new_cnv_site[1] = min(best_hit_cnv[1], pos)
                new_cnv_site[3] = max(best_hit_cnv[3], end)
                # ensure end does not exceed chrom length
                new_cnv_site[3] = min(new_cnv_site[3],
                                   chrom_length_dict[str(chrom)])
            # make DUP more specific, if possible
            if new_cnv_site[4] == "DUP":
                new_cnv_site[4] = sv_type
            # change new site to tuple
            new_cnv_site = tuple(new_cnv_site)
            # add new site to original sites
            new_cnv_original_sites.add(site)
            # create new cnv interval object
            new_site_dict = {new_cnv_site: set()}
            for original_site in new_cnv_original_sites:
                new_site_dict[new_cnv_site].add(original_site)
            # remove old cnv interval
            sites_interval_trees[chrom].remove(best_hit)
            # add new interval
            if new_cnv_site[4] == "INS":
                interval_start = max(0, pos - 11)
                interval_end = min(end + 10, chrom_length_dict[str(chrom)])
                sites_interval_trees[chrom].addi(interval_start, interval_end, new_site_dict)
            else:
                sites_interval_trees[chrom].addi(new_cnv_site[1] - 1, new_cnv_site[3], new_site_dict)
        else:
            # ensure that end does not exceed chromosome
            end = min(end, chrom_length_dict[str(chrom)])
            new_site = (chrom, pos, alt, end, sv_type, inschrom, inspos)
            new_site_dict = {new_site: set()}
            new_site_dict[new_site].add(new_site)
            # add new site
            if new_site[4] == "INS":
                # add 10 bp upstream and downstream for overlap queries
                interval_start = max(0, pos - 11)
                interval_end = min(end + 10, chrom_length_dict[str(chrom)])
                sites_interval_trees[chrom].addi(interval_start, interval_end, new_site_dict)
            else:
                sites_interval_trees[chrom].addi(new_site[1] - 1, new_site[3], new_site_dict)
    return sites_interval_trees

def write_to_output_vcf(output_vcf_fn, chrom_length_dict, sites_interval_trees, samples_dict):
    """
    Produce a VCF file containing all CNV sites and all genotypes of all samples

    :param output_vcf_fn: Name of VCF output file
    :param chrom_length_dict: Dictionary with chromosomes as keys and
    length of chromosomes as values
    :param sites_interval_tress: A dictionary containing chromosomes as keys
    and IntervalTree instances as values. The IntervalTree instances contain
    dictionaries with collapsed sites as keys and a set of non-collapsed sites
    as values. Sites are tuples
    (chrom, pos, alt, end, sv_type, inschrom, inspos).
    :param samples_dict: Dict containing sites as keys and genotype lines of
    samples as values.
    :return: 0 (integer)
    """
    # create identifier generator for sites
    sites_identifiers = count(1)
    sites_dict = {}
    for chrom in sites_interval_trees:
        for interval in sites_interval_trees[chrom]:
            sites_dict.update(interval.data)
    # sort sites, if not empty
    sorted_sites = []
    if sites_dict:
        sorted_sites = sorted(sites_dict.keys(), key=operator.itemgetter(0, 1, 3, 4, 2, 5, 6))
    # produce VCF header
    sample_names = samples_dict.keys()
    sample_names_header = "\t".join(sample_names)
    vcf_header_elems = ["##fileformat=VCFv4.3",
                        "##fileDate={}".format(datetime.date.today().strftime("%Y%m%d")),
                        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",
                        "##INFO=<ID=END,Number=.,Type=Integer,Description=\"End position of the variant described in this region\">",
                        "##INFO=<ID=INSCHROM,Number=.,Type=String,Description=\"Chromosome on which insertion site of the dispersed duplication is located\">",
                        "##INFO=<ID=INSPOS,Number=.,Type=Integer,Description=\"Position of insertion site of the dispersed duplication\">",
                        "##ALT=<ID=DEL,Description=\"Deletion\">",
                        "##ALT=<ID=DUP:DISPERSED,Description=\"Dispersed Duplication\">",
                        "##ALT=<ID=DUP:TANDEM,Description=\"Tandem Duplication\">",
                        "##ALT=<ID=INS,Description=\"Insertion of novel sequence\">",
                        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
                        "##FORMAT=<ID=SUP,Number=1,Type=Float,Description=\"Median number of supporting reads\">",
                        "##FORMAT=<ID=RP,Number=1,Type=Float,Description=\"Median number of supporting discordantly aligned read pairs\">",
                        "##FORMAT=<ID=SR,Number=1,Type=Float,Description=\"Median number of supporting split reads\">",
                        "##FORMAT=<ID=TOOL,Number=.,Type=String,Description=\"Supporting tools\">",
                        "##FORMAT=<ID=RQ,Number=1,Type=Float,Description=\"Probabilistic score of random forest model\">",
                        "##FORMAT=<ID=DHFC,Number=1,Type=Float,Description=\"duphold depth fold-change\">",
                        "##FORMAT=<ID=DHBFC,Number=1,Type=Float,Description=\"duphold depth fold-change compared to bins with matching GC\">",
                        "##FORMAT=<ID=DHFFC,Number=1,Type=Float,Description=\"duphold depth flank fold-change compared to 1000bp left and right of event\">"]
    for chrom in chrom_length_dict:
        contig_line = "##contig=<ID={},length={}>".format(chrom, chrom_length_dict[chrom])
        vcf_header_elems.append(contig_line)
    vcf_header_elems.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}".format(
        sample_names_header))
    vcf_header = "\n".join(vcf_header_elems)
    with open(output_vcf_fn, "w") as output_vcf:
        # write header
        output_vcf.write(vcf_header)
        output_vcf.write("\n")
        # loop through sites, writing information on each line for each site
        ref = "N"
        qual = "-1"
        filtered = "PASS"
        for collapsed_site in sorted_sites:
            chrom = str(collapsed_site[0])
            pos = str(collapsed_site[1])
            identifier = str(next(sites_identifiers))
            alt = collapsed_site[2]
            sv_type = collapsed_site[4]
            type_info_field = "SVTYPE={}".format(sv_type)
            info_field_elems = [type_info_field, ]
            if sv_type != "INS":
                end = "END={}".format(str(collapsed_site[3]))
                info_field_elems.append(end)
            # add insertion site for dispersed duplications
            if sv_type == "DUP:DISPERSED":
                ins_chrom = "INSCHROM={}".format(str(collapsed_site[5]))
                info_field_elems.append(ins_chrom)
                ins_pos = "INSPOS={}".format(str(collapsed_site[6]))
                info_field_elems.append(ins_pos)
            info_field = ";".join(info_field_elems)
            format_field = "GT:SUP:RP:SR:TOOL:RQ:DHFC:DHBFC:DHFFC"
            variant_line_elems = [chrom, pos, identifier, ref, alt, qual,
                                  filtered, info_field, format_field]
            # extract sample elements
            original_sites = sites_dict[collapsed_site]
            for sample_name in sample_names:
                sample_field = 0
                for original_site in original_sites:
                    if original_site in samples_dict[sample_name]:
                        sample_field = samples_dict[sample_name][original_site]
                if sample_field == 0:
                    sample_field_elems = ["0/0", "0", "0", "0", ".", "0", "-1", "-1", "-1"]
                    sample_field = ":".join(sample_field_elems)
                variant_line_elems.append(sample_field)
            # create new line for variant
            variant_line = "\t".join(variant_line_elems)
            output_vcf.write(variant_line)
            output_vcf.write("\n")
    return 0

def merge_vcfs(input_fn, output_vcf_fn, chrom_length_dict, reciprocal=0.5):
    """
    Merge VCF files of single samples containing CNVs called by Hecaton,
    producing a single VCF file containing the merged calls
    :param input_fn: Path to file containing path to single sample VCF file on each line
    :param output_vcf_fn: File containing path to VCF file
    :param chrom_length_dict: Dictionary with chromosomes as keys and
    length of chromosomes as values
    :param reciprocal: Minimum reciprocal overlap required to collapse CNVs
    :return: 0 (integer)
    """
    # get vcf filenames from input file
    input_fns = []
    with open(input_fn) as input_file:
        for line in input_file:
            input_fns.append(line.strip())
    # obtain all sites and genotype information
    sites_set, samples_dict = obtain_sites_and_genotypes(input_fns)
    # get interval tree in which sites considered to be the same CNV are collapsed
    sites_interval_tree = get_interval_tree_from_sites_set(sites_set,
                                                           chrom_length_dict,
                                                           reciprocal)
    # write sites and genotypes to output vcf
    write_to_output_vcf(output_vcf_fn, chrom_length_dict, sites_interval_tree, samples_dict)
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
    # check if reciprocal is correct
    if args.reciprocal < 0 or args.reciprocal > 1:
        raise ValueError("Reciprocal overlap must be between 0 and 1")
    # merge vcf files
    merge_vcfs(args.input_file, args.output_vcf, chrom_length_dict, args.reciprocal)

if __name__ == "__main__":
    main()

