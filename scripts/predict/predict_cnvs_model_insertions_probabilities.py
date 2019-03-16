#!/usr/bin/env python3

"""
Predict the probabilities of CNVs being true or false using a predictive model
"""

import pandas as pd
import pickle
import argparse
import sys

def parse_cl_args(in_args):
    """
    Parse command line arguments

    :param in_args: All command line arguments
    :return: None
    """
    description = "Predict the probability of a CNV being true or false " \
                  "from a set of extracted features"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-i", "--input_bedpe", type=str,
                        help="Path to bedpe file containing features")
    parser.add_argument("-m", "--model", type=str,
                        help="Path to pickle file containing model")
    parser.add_argument("-o", "--output_bedpe", type=str,
                        help="Name of bedpe output file")
    parser.add_argument("-s", "--scaler", type=str,
                        help="Name of scaler (optional)")
    args = parser.parse_args(in_args)
    return args

def predict_cnvs(input_bedpe, model_pickle_fn, output_bedpe, scaler_fn=None):
    """Predict cnvs from features

    :param input_bedpe: CSV file containing extracted features, produced by
    create_features_logistic_regression.py
    :param model_pickle_fn: Filehandle with reading rights to pickle file
    containing logistic model, produced by
    create_logistic_model_5x_solanum_lycopersicum.py
    :param output_bedpe: BEDPE file t0 which results will be written
    :param scaler_fn: Name of pickle file with scaler if working with
    a model that requires data to be scaled
    :return: 0
    """
    #load data frame
    input_df = pd.read_csv(input_bedpe, sep="\t")
    #load model
    model = pickle.load(model_pickle_fn)
    #extract features
    features = input_df.loc[:, ["READ_PAIRS", "SPLIT_READS", "SIZE", "DELLY", "LUMPY",
                                "MANTA", "GRIDSS", "DEL", "INS", "TANDUP", "DISDUP"]]
    # scale features if needed
    if scaler_fn:
        with open(scaler_fn, "rb") as scaler_file:
            scaler = pickle.load(scaler_file)
            features = scaler.transform(features)
    #print("Features: {0}".format(features))
    #predict cnvs, add nothing if bedpe is empty
    if not features.empty:
        preds = model.predict_proba(features)
        # add predictions to dataframe
        input_df["PREDICTION_0"] = preds[:, 0]
        input_df["PREDICTION_1"] = preds[:, 1]
    else:
        input_df["PREDICTION_0"] = ""
        input_df["PREDICTION_1"] = ""
    # write bedpe to output
    input_df.to_csv(output_bedpe, sep="\t", index=False)
    return 0

def main():
    args = parse_cl_args(sys.argv[1:])
    # predict cnvs
    with open(args.model, "rb") as model_file:
        predict_cnvs(args.input_bedpe, model_file, args.output_bedpe, args.scaler)

if __name__ == "__main__":
    main()