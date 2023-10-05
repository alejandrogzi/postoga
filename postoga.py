#!/usr/bin/env python3


__author__ = "Alejandro Gonzales-Irribarren"
__version__ = "0.1.0"
__email__ = "jose.gonzalesdezavala1@unmsm.edu.pe"
__github__ = "https://github.com/alejandrogzi"


import os
import pandas as pd
import numpy as np
import argparse
import sys
import subprocess
import time



def toga_table(args: argparse.Namespace) -> pd.DataFrame:
    """
    Return a pandas DataFrame with all projections and metadata.

    @type args: subprocess.Namespace
    @param args: defined arguments
    @rtype: DataFrame
    @return a DataFrame with all the projections and metada
    """
    # toga_table uses 3 sources within a TOGA directory: orthology_classification,
    # loss_sum_data, and ortholog_scores. These three files contain all the universe
    # of transcripts evaluated by TOGA. It initates by reading those three files as
    # pandas DataFrames.
    orthology = pd.read_csv(f"{args.path}/orthology_classification.tsv", sep="\t")
    loss = pd.read_csv(f"{args.path}/loss_summ_data.tsv", sep="\t", header=None, names=["projection", "transcript", "class"])
    score = pd.read_csv(f"{args.path}/temp/orthology_scores.tsv", sep="\t")
    isoforms = pd.read_csv(f"{args.path}/temp/isoforms.tsv", sep="\t", header=None)

    # Creates a dictionary: transcript -> gene
    isoforms_dict = isoforms.set_index(1).to_dict().get(0)
    
    # The loss DataFrame needs to be adjusted to only consider projections (TOGA predictions).
    # After filtering loss, toga_table() creates an additional helper column with the first part of 
    # each transcript.
    loss = loss[loss["projection"] == "PROJECTION"]
    loss["helper"] = loss["transcript"].str.split(".").str[0]

    ortho_x_loss = pd.merge(orthology, 
                            loss, 
                            left_on="q_transcript", 
                            right_on="transcript", 
                            how="outer")
    
    ortho_x_loss["helper"].fillna(ortho_x_loss["t_gene"], inplace=True)

    # Since the score DataFrame contains transcript names (under the column "gene") and chain
    # IDs (under the column "chain") in two separate columns, toga_table() joins both column to 
    # be able to merge scores in the next step
    score["transcripts"] = score["gene"] + "." + score["chain"].astype(str)
    score = score[["transcripts", "pred"]]

    # Merge ortho_x_loss and score DataFrames
    toga = pd.merge(ortho_x_loss, 
                    score, 
                    left_on="transcript", 
                    right_on="transcripts")

    # Get relevant column and fill missing gene names
    toga = toga[['t_gene', 'helper',  'transcripts', 'orthology_class','class', 'pred', 'q_gene']]
    toga["t_gene"].fillna(toga["helper"].map(isoforms_dict), inplace=True)
    
    return toga



def write_isoforms(args: argparse.Namespace, table: pd.DataFrame):
    """
    Writes all isoforms to a text file

    @type args: subprocess.Namespace
    @param args: defined arguments
    @type table: pd.DataFrame
    @param table: a pandas DataFrame
    """
    f = f"{args.path}/isoforms.txt"
    
    # Get only gene:transcript pairs
    table = table.iloc[:,[0,2]]
    table.to_csv(f, sep="\t", header=None, index=False)

    return



def filter_bed(args: argparse.Namespace, table: pd.DataFrame) -> str:
    """
    Filters the original .bed file to produce a custom filtered file

    @type args: subprocess.Namespace
    @param args: defined arguments
    @type table: pd.DataFrame
    @param table: a pandas DataFrame
    @rtype: str
    @return path: path to filtered file
    """

    if args.threshold:
        # Filter by threshold
        table = table[table["pred"] >= float(args.threshold)]

    if args.by_class:
        # Filter by class
        table = table[table["orthology_class"].isin(args.by_class.split(","))]

    if args.by_rel:
        # Filter by relationship
        table = table[table["class"].isin(args.by_rel.split(","))]


    # Read the original .bed file and filter it based on the transcripts table
    bed = pd.read_csv(f"{args.path}/query_annotation.bed", sep="\t", header=None)
    bed = bed[bed[3].isin(table["transcripts"])]

    # Write the filtered .bed file
    bed.to_csv(f"{args.path}/filtered.bed", sep="\t", header=None, index=False)

    return f"{args.path}/filtered.bed"



def run(args: argparse.Namespace):
    """
    The postoga runner function

    @type args: subprocess.Namespace
    @param args: defined arguments
    """
    # Get the toga table and write isoforms
    table = toga_table(args)
    write_isoforms(args, table)

    # Filtering step (if any)
    if any(args.by_class, args.by_rel, args.threshold):
        bed = filter_bed(args, table)
    else:
        bed = f"{args.path}/query_annotation.bed"




    return



def parser():
    app = argparse.ArgumentParser()
    app.add_argument(
        "-p",
        "--path",
        help="Path to TOGA results directory",
        required=True,
        type=str
    )
    app.add_argument(
        "-bc",
        "--by-class",
        help="Filter parameter to only include certain orthology classes (I, PI, UL, M, PM, L, UL)",
        required=False,
        type=str
    )
    app.add_argument(
        "-brel",
        "--by-rel",
        help="Filter parameter to only include certain orthology relationships (o2o, o2m, m2m, m2m, o2z)",
        required=False,
        type=str
    )
    app.add_argument(
        "-thold",
        "--threshold",
        help="Filter parameter to preserve orthology scores above a given treshold (0.0 - 1.0)",
        required=False,
        type=str
    )

    if len(sys.argv) < 1:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()

    return args


def main():
    args = parser()
    # run(args.path)


if __name__ == "__main__":
    main()