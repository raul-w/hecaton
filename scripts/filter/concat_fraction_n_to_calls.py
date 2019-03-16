#!/usr/bin/env python3

"""
Concatenate fraction of N's of flanking sequences to a BEDPE file with CNV calls
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
    description = "Add N fraction to detected BEDPE"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-i", "--detected_bedpe", type=str,
                        help="Path to bedpe file containing calls")
    parser.add_argument("-f", "--fraction_fn", type=str,
                        help="Path to file with fraction N's in breakend sequences")
    parser.add_argument("-o", "--output_bedpe", type=str,
                        help="Name of bedpe output file")
    args = parser.parse_args(in_args)
    return args

def concat_columns(input_bedpe, fraction_fn, output_bedpe):
    """Concatenate column with N fraction to training set

    :param input_bedpe: Path to bedpe file containing detected evaluations
    :param fraction_fn: Path to file with scores
    :param output_bedpe: Path to output bedpe file
    :return: 0
    """
    # load data frames
    input_df = pd.read_csv(input_bedpe, sep="\t")
    fraction_df = pd.read_csv(fraction_fn, sep="\t")
    if input_df.empty:
        # add empty column
        input_df["Flanking_Ns"] = ""
        # write dataframe to bedpe file
        input_df.to_csv(output_bedpe, sep="\t", index=False)
        return 0
    # add alignment score to dataframe
    input_df["Flanking_Ns"] = fraction_df["Fraction_Ns"]
    # write bedpe to output
    input_df.to_csv(output_bedpe, sep="\t", index=False)
    return 0

def main():
    args = parse_cl_args(sys.argv[1:])
    # concatenate column
    concat_columns(args.detected_bedpe, args.fraction_fn, args.output_bedpe)

if __name__ == "__main__":
    main()