#!/usr/bin/env python3

"""
Calculate fraction of N's in breakend sequences identified by bedtools
"""

import argparse
import sys
from Bio import Align, SeqIO

def parse_cl_args(in_args):
    """
    Parse command line arguments

    :param in_args: All command line arguments
    :return: None
    """
    description = "Calculate fraction of N's in breakend sequences identified by bedtools"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-i", "--breakend_fasta", type=str,
                        help="FASTA file containing breakend sequences")
    parser.add_argument("-o", "--output_fn", type=str,
                        help="Name of output file")
    args = parser.parse_args(in_args)
    return args

def fraction_ns(seq1: str, seq2: str):
    """Calculate fraction of N's in a pair of DNA sequences

    :param seq1: DNA sequence (str)
    :param seq2: DNA sequence (str)
    :return: Fraction of sequences consisting of N's (float)
    """
    # change sequence to upper case
    seq1 = seq1.upper()
    seq2 = seq2.upper()
    # check if both are DNA strings
    allowed_bases = "ATCGNKMRYSWBVHDX.-"
    for ch in seq1:
        if ch not in allowed_bases:
            raise ValueError("Not a DNA FASTA sequence: contains {}".format(ch))
    for ch in seq2:
        if ch not in allowed_bases:
            raise ValueError("Not a DNA FASTA sequence: contains {}".format(ch))
    # calculate fraction of N's
    seq1_ns = seq1.count("N")
    seq2_ns = seq2.count("N")
    fraction_n = (seq1_ns + seq2_ns) / (len(seq1) + len(seq2))
    return fraction_n

def fraction_ns_breakend_file(breakend_fasta_fn, output_fn):
    """Calculate fraction of N's in breakend sequences identified by bedtools

    :param breakend_fasta_fn: Breakend FASTA file
    :param output_fn: Output file
    :return: 0 (integer)
    """
    # initialize list of N fractions
    fractions = []
    # calculate fraction of N's for all pairs of breakend sequences
    # and store the fraction
    with open(breakend_fasta_fn) as multiFASTA:
        records = SeqIO.parse(multiFASTA, 'fasta')
        fraction_time = False
        seq1 = ""
        seq2 = ""
        for record in records:
            if not fraction_time:
                seq1 = record.seq
                fraction_time = True
            elif fraction_time:
                seq2 = record.seq
                fraction = fraction_ns(seq1, seq2)
                fractions.append(fraction)
                fraction_time = False
    # write scores to output file
    with open(output_fn, "w") as output_file:
        output_file.write("Fraction_Ns\n")
        score_template = "{}\n"
        for fraction in fractions:
            output_file.write(score_template.format(fraction))
    return 0

def main():
    args = parse_cl_args(sys.argv[1:])
    # get fractions of N's
    fraction_ns_breakend_file(args.breakend_fasta, args.output_fn)

if __name__ == "__main__":
    main()