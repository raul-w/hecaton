#!/usr/bin/env python3

"""
Convert BEDPE files containing intersected sets of CNVs from different tools
to a BEDPE file containing the features used for a logistic model
"""

import argparse
import pandas as pd
import sys
from statistics import median

def parse_cl_args(in_args):
    """
    Parse command line arguments

    :param in_args: All command line arguments
    :return: None
    """
    description = "Convert intersect BEDPE to feature BEDPE"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-i", "--input_bedpe", type=str,
                        help="Path to intersect BEDPE file")
    parser.add_argument("-o", "--output_bedpe", type=str,
                        help="Name of feature BEDPE file")
    args = parser.parse_args(in_args)
    return args

def intersect_to_feature(input_bedpe, output_bedpe):
    """Convert intersect BEDPE to feature BEDPE

    :param input_bedpe: Filename of intersect BEDPE file
    :param output_bedpe: Filename of feature BEDPE file
    :param interaction: Boolean specifying if interaction features should be added
    :return: 0 (integer)
    """
    # load input bedpe
    input_df = pd.read_csv(input_bedpe, sep="\t")
    # check if input_df is empty
    if input_df.empty:
        # add empty columns
        new_column_headers = ["SIZE", "READ_PAIRS", "SPLIT_READS", "CNVNATOR",
                              "CONTROL-FREEC", "DELLY", "LUMPY", "MANTA", "GRIDSS",
                              "DEL", "INS", "DUP", "TANDUP", "DISDUP"]
        # remove original tool column
        input_df.drop("TOOL", axis=1, inplace=True)
        input_df = input_df.reindex(columns=input_df.columns.tolist() + new_column_headers)
        # write dataframe to bedpe file
        input_df.to_csv(output_bedpe, sep="\t", index=False)
        return 0
    # create size column
    size_column = []
    # create read pairs and split reads columns
    read_pair_column = []
    split_read_column = []
    # create tool columns
    cnvnator_column = []
    freec_column = []
    delly_column = []
    lumpy_column = []
    manta_column = []
    gridss_column = []
    # create sv type columns
    del_column = []
    ins_column = []
    dup_column = []
    tandup_column = []
    disdup_column = []
    # check which cnvs were detected by which tools
    for predicted_cnv in input_df.itertuples(index=False, name='Pandas'):
        # initialize tool prediction booleans
        cnvnator_detected = 0
        freec_detected = 0
        delly_detected = 0
        lumpy_detected = 0
        manta_detected = 0
        gridss_detected = 0
        tools = predicted_cnv.TOOL.split(";")
        for tool in tools:
            tool = tool.upper()
            if tool == "CNVNATOR":
                cnvnator_detected = 1
            if tool == "CONTROL-FREEC":
                freec_detected = 1
            if tool == "DELLY":
                delly_detected = 1
            if tool == "LUMPY":
                lumpy_detected = 1
            if tool == "MANTA":
                manta_detected = 1
            if tool == "GRIDSS":
                gridss_detected = 1
        # append tool booleans to columns
        cnvnator_column.append(cnvnator_detected)
        freec_column.append(freec_detected)
        delly_column.append(delly_detected)
        lumpy_column.append(lumpy_detected)
        manta_column.append(manta_detected)
        gridss_column.append(gridss_detected)
        # initialize sv type booleans
        deletion = 0
        insertion = 0
        dup = 0
        tandup = 0
        disdup = 0
        # fill in sv type
        sv_type = predicted_cnv.TYPE
        if sv_type == "DEL":
            deletion = 1
        elif sv_type == "INS":
            insertion = 1
        elif sv_type == "DUP":
            dup = 1
        elif sv_type == "DUP:TANDEM":
            tandup = 1
        elif sv_type == "DUP:DISPERSED":
            disdup = 1
        del_column.append(deletion)
        ins_column.append(insertion)
        dup_column.append(dup)
        tandup_column.append(tandup)
        disdup_column.append(disdup)
        # create correct size
        if sv_type == "INS":
            insertion_sequence = predicted_cnv.INSERTED_SEQUENCE
            if insertion_sequence == "N":
                size = 5000
            else:
                size = len(insertion_sequence)
        else:
            size = predicted_cnv.END_A - predicted_cnv.START_A
        size_column.append(size)
        # take median of read pairs and split reads
        read_pairs = [float(i) for i in predicted_cnv.READ_PAIRS.split(";")]
        read_pair_column.append(median(read_pairs))
        split_reads = [float(i) for i in predicted_cnv.SPLIT_READS.split(";")]
        split_read_column.append(median(split_reads))
    # add size column to df
    input_df["SIZE"] = size_column
    # replace read pairs and split reads columns
    input_df["READ_PAIRS"] = read_pair_column
    input_df["SPLIT_READS"] = split_read_column
    # add tool columns to df
    input_df["CNVNATOR"] = cnvnator_column
    input_df["CONTROL-FREEC"] = freec_column
    input_df["DELLY"] = delly_column
    input_df["LUMPY"] = lumpy_column
    input_df["MANTA"] = manta_column
    input_df["GRIDSS"] = gridss_column
    # remove original tool column
    input_df.drop("TOOL", axis=1, inplace=True)
    # add sv type columns to df
    input_df["DEL"] = del_column
    input_df["INS"] = ins_column
    input_df["DUP"] = dup_column
    input_df["TANDUP"] = tandup_column
    input_df["DISDUP"] = disdup_column
    # write dataframe to bedpe file
    input_df.to_csv(output_bedpe, sep="\t", index=False)
    return 0

def main():
    args = parse_cl_args(sys.argv[1:])
    # convert vcf file
    intersect_to_feature(args.input_bedpe, args.output_bedpe)

if __name__ == "__main__":
    main()