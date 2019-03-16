#!/usr/bin/env python3

"""
Fix the VCF header of GRIDSS files
"""

import argparse
from xopen import xopen
import sys

def parse_cl_args(in_args):
    """
    Parse command line arguments

    :param in_args: All command line arguments
    :return: None
    """
    description = "Fix the VCF header of GRIDSS files"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-i", "--input_vcf_fn", type=str,
                        help="Path to VCF file")
    parser.add_argument("-o", "--output_vcf_fn", type=str,
                        help="Name of output VCF file")
    args = parser.parse_args(in_args)
    return args

def fix_vcf_header(input_vcf, output_vcf):
    """Fix the VCF header of GRIDSS files

    :param input_vcf: Filehandle with reading rights to GRIDSS VCF file
    :param output_vcf: Filehandle with writing rights to GRIDSS VCF file
    :return: 0 (integer)
    """
    # fix header
    for line in input_vcf:
        if line.startswith("##"):
            if "BANRPQ" in line or "BANSRQ" in line:
                line = line.replace("Integer", "Float", 1)
            output_vcf.write(line)
        else:
            output_vcf.write(line)
    return 0

def main():
    args = parse_cl_args(sys.argv[1:])
    # fix header
    with xopen(args.input_vcf_fn) as input_vcf, \
        xopen(args.output_vcf_fn, "w") as output_vcf:
        fix_vcf_header(input_vcf, output_vcf)

if __name__ == "__main__":
    main()