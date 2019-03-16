#!/usr/bin/env python3

"""
Convert BEDPE file containing CNV calls to VCF
"""

import argparse
import datetime
import pandas as pd
import sys


def parse_cl_args(in_args):
    """
    Parse command line arguments

    :param in_args: All command line arguments
    :return: None
    """
    description = "Convert BEDPE file containing CNV calls to VCF"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-i", "--input_bedpe", type=str,
                        help="Path to BEDPE file")
    parser.add_argument("-o", "--output_vcf", type=str,
                        help="Name of VCF output file")
    parser.add_argument("-s", "--sample_name", type=str,
                        help="Name of the sample")
    args = parser.parse_args(in_args)
    return args


def convert_bedpe(input_bedpe_fn, output_vcf_fn, sample_name):
    """Convert BEDPE file containing CNV calls to VCF

    :param input_bedpe_fn: Path to input BEDPE file
    :param output_vcf_fn: Path to output BEDPE file
    :param sample_name: Name of the sample
    :return: 0 (integer)
    """
    # load bedpe file
    input_bedpe = pd.read_csv(input_bedpe_fn, sep="\t")
    # write VCF header to output file
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
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}".format(sample_name)]
    vcf_header = "\n".join(vcf_header_elems)
    with open(output_vcf_fn, "w") as output_vcf:
        output_vcf.write(vcf_header)
        output_vcf.write("\n")
        # loop through bedpe file writing each variant to the vcf
        for cnv in input_bedpe.itertuples(index=False, name='Pandas'):
            chrom = str(cnv[0])
            sv_type = cnv.TYPE
            if sv_type == "INS":
                start = str(cnv.START_A)
            else:
                start = str(cnv.START_A + 1)  # vcf format
            identifier = str(cnv.ID)
            ref = "N"
            if sv_type == "INS":
                alt = "".join(["N", cnv.INSERTED_SEQUENCE])
            else:
                alt = "<{}>".format(sv_type.split(":")[0])
            qual = str(cnv.QUAL)
            filtered = "PASS"
            type_info_field = "SVTYPE={}".format(sv_type)
            info_field_elems = [type_info_field, ]
            if sv_type != "INS":
                end = "END={}".format(str(cnv.END_A))
                info_field_elems.append(end)
            # add insertion site for dispersed duplications
            if sv_type == "DUP:DISPERSED":
                ins_chrom = "INSCHROM={}".format(str(cnv.CHROM_B))
                info_field_elems.append(ins_chrom)
                ins_pos = "INSPOS={}".format(cnv.START_B)
                info_field_elems.append(ins_pos)
            info_field = ";".join(info_field_elems)
            # extract format field elements
            format_field = "GT:SUP:RP:SR:TOOL:RQ"
            read_pairs = float(cnv.READ_PAIRS)
            split_reads = float(cnv.SPLIT_READS)
            total_support = read_pairs + split_reads
            # extract tools
            tools = []
            if cnv.DELLY:
                tools.append("DELLY")
            if cnv.LUMPY:
                tools.append("LUMPY")
            if cnv.MANTA:
                tools.append("MANTA")
            if cnv.GRIDSS:
                tools.append("GRIDSS")
            tool_field = ",".join(tools)
            rf_score = round(cnv.PREDICTION_1, 2)
            sample_field_elems = ["1/1", str(total_support), str(read_pairs),
                                  str(split_reads), tool_field, str(rf_score)]
            sample_field = ":".join(sample_field_elems)
            # create new line for variant
            variant_line_elems = [chrom, start, identifier, ref, alt, qual,
                                  filtered, info_field, format_field, sample_field]
            variant_line = "\t".join(variant_line_elems)
            output_vcf.write(variant_line)
            output_vcf.write("\n")
    return 0


def main():
    args = parse_cl_args(sys.argv[1:])
    # convert vcf file
    convert_bedpe(args.input_bedpe, args.output_vcf, args.sample_name)


if __name__ == "__main__":
    main()
