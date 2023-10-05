#!/usr/bin/env python3

"""
Master script for postoga. 

postoga is a tool that automates the post-processing of TOGA results.
At its core, this tool takes a TOGA results directory and produces
a series of steps to reduce the amount of manual work required to pre-process
files for downstream analysis.
"""


__author__ = "Alejandro Gonzales-Irribarren"
__version__ = "0.1.0-devel"
__email__ = "jose.gonzalesdezavala1@unmsm.edu.pe"
__github__ = "https://github.com/alejandrogzi"



import os
import pandas as pd
import numpy as np
import argparse
import sys
import subprocess
import time
from modules.constants import Constants
# from version import __version__



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

    orthology = pd.read_csv(os.path.join(args.path, Constants.FileNames.ORTHOLOGY), sep="\t")
    loss = pd.read_csv(os.path.join(args.path, Constants.FileNames.CLASS), sep="\t", header=None, names=["projection", "transcript", "class"])
    score = pd.read_csv(os.path.join(args.path, Constants.FileNames.SCORES), sep="\t")
    isoforms = pd.read_csv(os.path.join(args.path, Constants.FileNames.ISOFORMS), sep="\t", header=None)

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
    f = os.path.join(args.path, Constants.FileNames.OWNED_ISOFORMS)
    
    # Get only gene:transcript pairs
    table = table.iloc[:,[0,2]]
    table.to_csv(f, sep="\t", header=None, index=False)

    return f



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
    bed = pd.read_csv(os.path.join(args.path, Constants.FileNames.BED), sep="\t", header=None)
    bed = bed[bed[3].isin(table["transcripts"])]

    # Write the filtered .bed file
    bed.to_csv(os.path.join(args.path, Constants.FileNames.FILTERED_BED), sep="\t", header=None, index=False)

    return os.path.join(args.path, Constants.FileNames.FILTERED_BED)



def shell(cmd):
    """ 
    Run a shell command and return the output as a string
    
    @type cmd: str
    @param cmd: shell command
    """
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return result.stdout.strip()



def bed_to_gtf(bed: str, isoforms: str):
    """
    Converts a .bed file to .gtf

    @type bed: str
    @param bed: path to .bed file
    @type isoforms: str
    @param isoforms: path to the isoforms file
    """

    gtf = f"{bed.split('.')[0]}.gtf"
    cmd = f"{Constants.ToolNames.BED2GTF} {bed} {isoforms} {gtf}"
    sh = shell(cmd)

    return gtf



def bed_to_gff(bed: str, isoforms: str):
    """
    Converts a .bed file to .gff

    @type bed: str
    @param bed: path to .bed file
    @type isoforms: str
    @param isoforms: path to the isoforms file
    """
    gff = f"{bed.split('.')[0]}.gff"
    cmd = f"{Constants.ToolNames.BED2GFF} {bed} {isoforms} {gff}"
    sh = shell(cmd)

    return gff



def run(args: argparse.Namespace):
    """
    The postoga runner function

    @type args: subprocess.Namespace
    @param args: defined arguments
    """

    # Get the toga table and write isoforms
    table = toga_table(args)
    isoforms = write_isoforms(args, table)

    # Filtering step (if any)
    if any([args.by_class, args.by_rel, args.threshold]):
        bed = filter_bed(args, table)
    else:
        bed = f"{args.path}/query_annotation.bed"

    # Conversion step
    if args.to == "gtf":
        bed_to_gtf(bed, isoforms)
    elif args.to == "gff":
        bed_to_gff(bed, isoforms)


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
        help="Filter parameter to preserve orthology scores greater or equal a given threshold (0.0 - 1.0)",
        required=False,
        type=str
    )
    app.add_argument(
        "-to",
        "--to",
        help="Specify the conversion format for .bed (query_annotation/filtered) file (gtf, gff3)",
        required=True,
        type=str,
        choices=["gtf", "gff"]
    )

    if len(sys.argv) < 1:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()

    return args


def main():
    args = parser()
    run(args)


if __name__ == "__main__":
    main()