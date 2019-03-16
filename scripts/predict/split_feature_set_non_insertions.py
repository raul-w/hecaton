#!/usr/bin/env python3

"""
Split feature set into insertions and non-insertions
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
    description = "Split feature set into insertions and non-insertions"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-f", "--feature_bedpe", type=str,
                        help="Path to bedpe file containing features")
    parser.add_argument("-i", "--output_insertion_bedpe", type=str,
                        help="Name of output file containing insertions")
    parser.add_argument("-n", "--output_non_insertion_bedpe", type=str,
                        help="Name of output file containing non-insertions")
    args = parser.parse_args(in_args)
    return args

def filter_insertions(feature_bedpe, insertion_bedpe, non_insertion_bedpe):
    """Split feature bedpe into a bedpe containing insertions and a bedpe
    containing the rest

    :param feature_bedpe: Path to bedpe file containing calls
    :param insertion_bedpe: Path to output bedpe file containing insertions
    :param non_insertion_bedpe: Path to output bedpe file containing
    non-insertions
    :return: 0
    """
    # load data frame
    input_df = pd.read_csv(feature_bedpe, sep="\t")
    if input_df.empty:
        insertion_df = input_df
        non_insertion_df = input_df
    else:
        # split set into insertions and non-insertions
        insertion_df = input_df[input_df["TYPE"] == "INS"]
        non_insertion_df = input_df[input_df["TYPE"] != "INS"]
    # write bedpe to output
    insertion_df.to_csv(insertion_bedpe, sep="\t", index=False)
    non_insertion_df.to_csv(non_insertion_bedpe, sep="\t", index=False)
    return 0

def main():
    args = parse_cl_args(sys.argv[1:])
    # concatenate column
    filter_insertions(args.feature_bedpe, args.output_insertion_bedpe,
                      args.output_non_insertion_bedpe)

if __name__ == "__main__":
    main()