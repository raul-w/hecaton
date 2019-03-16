#!/usr/bin/env python3

"""
Convert VCF files to site-only VCF, removing all genotype information
"""

import argparse
import sys
from xopen import xopen


def parse_cl_args(in_args):
    """
    Parse command line arguments

    :param in_args: All command line arguments
    :return: None
    """
    description = "Convert VCF file to tabular file"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-v", "--vcf_fn", type=str,
                        help="File containing path to VCF file")
    parser.add_argument("-o", "--output_fn", type=str,
                        help="Name of VCF output file")
    args = parser.parse_args(in_args)
    return args

def vcf_to_site_only_vcf(input_fn, output_fn):
    """
    Convert VCF files to site-only VCF, removing all genotype information

    :param input_fn: Path to VCF file
    :param output_fn: Name of VCF output file
    :return: 0 (integer)
    """
    with xopen(input_fn) as input_file, xopen(output_fn, "w") as output_file:
        for line in input_file:
            # write initial header lines to output file
            if line.startswith("##"):
                output_file.write(line)
            else:
                # strip samples from line
                output_line = '\t'.join(line.strip().split()[0:9])
                output_file.write(output_line)
                output_file.write("\n")
    return 0

def main():
    args = parse_cl_args(sys.argv[1:])
    # write VCF fields to tabular output
    vcf_to_site_only_vcf(args.vcf_fn, args.output_fn)

if __name__ == "__main__":
    main()
