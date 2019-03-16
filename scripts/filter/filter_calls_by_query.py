#!/usr/bin/env python3

"""
Filter CNVs from a BEDPE file using a a query with bitwise operators
"""

import pandas as pd
import argparse
import sys

def parse_cl_args(in_args):
    """
    Parse command line arguments

    :param in_args: All command line arguments
    :return: None
    """
    description = "Add VaPoR score to detected BEDPE"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-i", "--detected_bedpe", type=str,
                        help="Path to bedpe file containing detected evaluations")
    parser.add_argument("-q", "--query", type=str,
                        help="Query to use for filtering")
    parser.add_argument("-o", "--output_bedpe", type=str,
                        help="Name of bedpe output file")
    args = parser.parse_args(in_args)
    return args

def filter_rows(input_bedpe, query, output_bedpe):
    """Filter CNVs from a BEDPE file using a query with bitwise operators

    :param input_bedpe: Path to bedpe file containing calls
    :param score: Score threshold
    :param output_bedpe: Path to output bedpe file
    :return: 0
    """
    # load data frame
    input_df = pd.read_csv(input_bedpe, sep="\t")
    # keep calls which have a score above the threshold
    input_df = input_df.query(query)
    # write bedpe to output
    input_df.to_csv(output_bedpe, sep="\t", index=False)
    return 0

def main():
    args = parse_cl_args(sys.argv[1:])
    # concatenate column
    filter_rows(args.detected_bedpe, args.query, args.output_bedpe)

if __name__ == "__main__":
    main()