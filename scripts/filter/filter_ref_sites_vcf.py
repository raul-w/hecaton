#!/usr/bin/env python3

"""
Filter sites in a VCF file for which the number of samples carrying a
variant allele goes beyond a user-chosen threshold
"""

import argparse
import pysam
import sys


def parse_cl_args(in_args):
    """
    Parse command line arguments

    :param in_args: All command line arguments
    :return: None
    """
    description = "Filter sites in a VCF file for which a minimum number of samples were found carrying a variant allele"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-v", "--vcf_fn", type=str,
                        help="File containing path to VCF file")
    parser.add_argument("-n", "--minimum_samples", type=int,
                        help="Minimum number of samples that should carry a variant allele")
    parser.add_argument("-o", "--output_fn", type=str,
                        help="Name of output file")
    args = parser.parse_args(in_args)
    return args

def filter_ref_sites(input_fn, min_samples, output_fn):
    """
    Filter sites in a VCF file for which the number of samples carrying a
    variant allele goes beyond a user-chosen threshold

    :param input_fn: Path to VCF file
    :param min_samples: Minimum number of samples carrying a variant allele
    :param output_fn: Name of VCF output file
    :return: 0 (integer)
    """
    with pysam.VariantFile(input_fn) as vcf, open(output_fn, 'w') as output_file:
        # get sample names
        samples = list(vcf.header.samples)
        # write header to output file
        header = str(vcf.header)
        output_file.write(header)
        # only write records which have non-reference calls
        for record in vcf.fetch():
            var_calls = 0
            # compute total number of variant calls
            if "GT" not in record.format:
                raise ValueError("GT not in format")
            else:
                for sample in samples:
                    genotype = record.samples[sample]["GT"]
                    print(genotype)
                    non_variants = [(0, 0), (None, None), (None, 0)]
                    if genotype not in non_variants:
                        var_calls += 1
            # write record to output if it has variant calls
            if var_calls >= min_samples:
                output_file.write(str(record))
    return 0

def main():
    args = parse_cl_args(sys.argv[1:])
    # filter ref sites from input file
    filter_ref_sites(args.vcf_fn, args.minimum_samples, args.output_fn)

if __name__ == "__main__":
    main()
