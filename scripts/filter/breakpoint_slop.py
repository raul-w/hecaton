#!/usr/bin/env python3

"""
Add slop to breakpoint sequences of CNVs, taking into account small CNVs with
a size smaller than the slop
"""

import argparse
import sys
import pandas as pd

def parse_cl_args(in_args):
    """
    Parse command line arguments

    :param in_args: All command line arguments
    :return: None
    """
    description = "Add slop to breakpoint sequences of CNVs, taking into " \
                  "account small CNVs with a size smaller than the slop"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-i", "--input_fn", type=str,
                        help="Name of input BEDPE file")
    parser.add_argument("-o", "--output_fn", type=str,
                        help="Name of output BEDPE file")
    parser.add_argument("-s", "--slop", type=int,
                        help="Amount of bp slop to add")
    args = parser.parse_args(in_args)
    return args

def compute_slop(ref_size, max_slop_size):
    """
    Compute slop to add to breakpoint sequences of CNVs, taking into
    account size of reference sequence affected by CNV

    :param ref_size: Size of sequence in reference affected by CNV
    :param max_slop_size: Maximum size of slop
    :return: Amount of slop (integer)
    """
    # raise error if max_slop_size or ref_size is negative
    if ref_size < 0 or max_slop_size < 0:
        raise ValueError("Ref size or max slop size is lower than 0")
    # compute maximum amount of slop that can be added without introducing
    # overlap
    slop_limit = int(ref_size / 2)
    slop = min(slop_limit, max_slop_size)
    return slop

def add_slop(input_fn, output_fn, max_slop):
    """Add slop to breakpoint sequences of CNVs, taking into
    account small CNVs with a size smaller than the slop

    :param input_fn: Name of input BEDPE file
    :param output_fn: Name of output BEDPE file
    :param max_slop: Maximum size of slop that needs to be added (int)
    :return: 0 (integer)
    """
    # load data frame
    input_df = pd.read_csv(input_fn, sep="\t")
    # get size of reference sequence affected by CNV
    ends = input_df["END_A"]
    starts = input_df["START_A"]
    ref_sizes = ends - starts
    # add slop to starts and ends
    new_starts = []
    new_ends = []
    for i, size in enumerate(ref_sizes):
        slop = compute_slop(size, max_slop)
        new_starts.append(starts[i] + slop)
        new_ends.append(ends[i] - slop)
    # replace starts and ends
    input_df["START_A"] = new_starts
    input_df["END_A"] = new_ends
    # write to new output file
    input_df.to_csv(output_fn, sep="\t", index=False)
    return 0

def main():
    args = parse_cl_args(sys.argv[1:])
    # add slop
    add_slop(args.input_fn, args.output_fn, args.slop)

if __name__ == "__main__":
    main()