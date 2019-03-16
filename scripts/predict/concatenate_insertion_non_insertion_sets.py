#!/usr/bin/env python

"""
Concatenate insertion and non-insertion sets
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
    description = "Concatenate insertion and non-insertion sets"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-i", "--insertion_bedpe", type=str,
                        help="Path to insertion set")
    parser.add_argument("-n", "--non_insertion_bedpe", type=str,
                        help="Path to non-insertion set")
    parser.add_argument("-o", "--output_bedpe", type=str,
                        help="Name of bedpe output file")
    args = parser.parse_args(in_args)
    return args

def concat_sets(insertion_fn, non_insertion_fn, output_fn):
    """Concatenate insertion and non-insertion sets

    :param insertion_fn: Filename of insertion training set
    :param non_insertion_fn: Filename of non-insertion training set
    :param output_fn: Filename of output file
    :return: 0 (integer)
    """
    # load the dataframes
    dataframes = []
    insertion_df = pd.read_csv(insertion_fn, sep="\t")
    dataframes.append(insertion_df)
    non_insertion_df = pd.read_csv(non_insertion_fn, sep="\t")
    dataframes.append(non_insertion_df)
    # concatenate the dataframes
    res = pd.concat(dataframes, ignore_index=True)
    # sort the dataframe
    res.sort_values(by=["#CHROM_A", "START_A", "END_A"], inplace=True)
    res.to_csv(output_fn, sep="\t", index=False)
    return 0

def main():
    args = parse_cl_args(sys.argv[1:])
    # concatenate column
    concat_sets(args.insertion_bedpe, args.non_insertion_bedpe, args.output_bedpe)

if __name__ == "__main__":
    main()