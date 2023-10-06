#!/usr/bin/env python3


"""
Master script for postoga. 

postoga is a tool that automates the post-processing of TOGA results.
At its core, this tool takes a TOGA results directory and produces
a series of steps to reduce the amount of manual work required to pre-process
files for downstream analysis.
"""


import os
import pandas as pd
import numpy as np
import argparse
import sys
import subprocess
import time
from constants import Constants
from version import __version__
from logger import Log


__author__ = "Alejandro Gonzales-Irribarren"
__version__ = __version__
__email__ = "jose.gonzalesdezavala1@unmsm.edu.pe"
__github__ = "https://github.com/alejandrogzi"


log = None


def toga_table(args: argparse.Namespace) -> pd.DataFrame:
    """
    Return a pandas DataFrame with all projections and metadata.

    @type args: subprocess.Namespace
    @param args: defined arguments
    @rtype: DataFrame
    @return a DataFrame with all the projections and metada
    """

    global log

    # Reads orthology_classification, loss_sum_data, and ortholog_scores.
    orthology = pd.read_csv(
        os.path.join(args.path, Constants.FileNames.ORTHOLOGY), sep="\t"
    )
    loss = pd.read_csv(
        os.path.join(args.path, Constants.FileNames.CLASS),
        sep="\t",
        header=None,
        names=["projection", "transcript", "class"],
    )
    score = pd.read_csv(os.path.join(args.path, Constants.FileNames.SCORES), sep="\t")
    isoforms = pd.read_csv(
        os.path.join(args.path, Constants.FileNames.ISOFORMS), sep="\t", header=None
    )
    quality = pd.read_csv(
        os.path.join(args.path, Constants.FileNames.QUALITY), sep="\t"
    )

    # Creates a dictionary: transcript -> gene
    isoforms_dict = isoforms.set_index(1).to_dict().get(0)

    # Subsets loss to consider only projections
    loss = loss[loss["projection"] == "PROJECTION"]
    loss["helper"] = loss["transcript"].str.split(".").str[0]

    ortho_x_loss = pd.merge(
        orthology, loss, left_on="q_transcript", right_on="transcript", how="outer"
    )

    ortho_x_loss["helper"].fillna(ortho_x_loss["t_gene"], inplace=True)

    # Merge transcript names (under the column "gene") and chain IDs (under the column "chain")
    score["transcripts"] = score["gene"] + "." + score["chain"].astype(str)
    score = score[["transcripts", "pred"]]

    table = pd.merge(ortho_x_loss, score, left_on="transcript", right_on="transcripts")

    # Create a new column with a rename orthology relationship
    table["relation"] = table["orthology_class"].map(Constants.ORTHOLOGY_TYPE)
    table["t_gene"].fillna(table["helper"].map(isoforms_dict), inplace=True)

    # Merge quality data
    table = pd.merge(table, quality, left_on="transcripts", right_on="Projection_ID")

    table = table[
        [
            "t_gene",
            "helper",
            "transcripts",
            "relation",
            "class",
            "pred",
            "q_gene",
            "confidence_level",
        ]
    ]

    log.record(
        f"found {len(table)} projections, {len(table['helper'].unique())} unique transcripts, {len(table['t_gene'].unique())} unique genes"
    )
    log.record(f"class stats: {table['class'].value_counts().to_dict()}")
    log.record(f"relation stats: {table['relation'].value_counts().to_dict()}")
    log.record(
        f"confidence stats: {table['confidence_level'].value_counts().to_dict()}"
    )

    return table


def write_isoforms(args: argparse.Namespace, table: pd.DataFrame) -> str:
    """
    Writes all isoforms to a text file

    @type args: subprocess.Namespace
    @param args: defined arguments
    @type table: pd.DataFrame
    @param table: a pandas DataFrame
    """

    global log

    f = os.path.join(args.path, Constants.FileNames.OWNED_ISOFORMS)

    # Get only gene:transcript pairs
    table = table.iloc[:, [0, 2]]
    table.to_csv(f, sep="\t", header=None, index=False)

    log.record(f"gene-to-projection hash with {len(table)} entries written to {f}")

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

    global log

    initial = len(table)

    if args.threshold:
        table = table[table["pred"] >= float(args.threshold)]
        log.record(
            f"discarded {initial - len(table)} projections with orthology scores <{args.threshold}"
        )

    if args.by_class:
        edge = len(table)
        table = table[table["class"].isin(args.by_class.split(","))]
        log.record(
            f"discarded {edge - len(table)} projections with classes other than {args.by_class}"
        )

    if args.by_rel:
        edge = len(table)
        table = table[table["relation"].isin(args.by_rel.split(","))]
        log.record(
            f"discarded {edge - len(table)} projections with relationships other than {args.by_rel}"
        )

    # Read the original .bed file and filter it based on the transcripts table
    bed = pd.read_csv(
        os.path.join(args.path, Constants.FileNames.BED), sep="\t", header=None
    )
    bed = bed[bed[3].isin(table["transcripts"])]
    custom_table = table[table["transcripts"].isin(bed[3])]

    # Write the filtered .bed file
    f = os.path.join(args.path, Constants.FileNames.FILTERED_BED)
    bed.to_csv(
        f,
        sep="\t",
        header=None,
        index=False,
    )

    log.record(
        f"kept {len(bed)} projections after filters, discarded {initial - len(bed)}."
    )
    log.record(
        f"{len(bed)} projections are coming from {len(custom_table['helper'].unique())} unique transcripts and {len(custom_table['t_gene'].unique())} genes"
    )
    log.record(
        f"class stats of new bed: {custom_table['class'].value_counts().to_dict()}"
    )
    log.record(
        f"relation stats of new bed: {custom_table['relation'].value_counts().to_dict()}"
    )
    log.record(
        f"confidence stats of new bed: {custom_table['confidence_level'].value_counts().to_dict()}"
    )
    log.record(f"filtered bed file written to {f}")

    return f


def shell(cmd) -> str:
    """
    Run a shell command and return the output as a string

    @type cmd: str
    @param cmd: shell command
    """
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return result.stdout.strip()


def bed_to_gtf(bed: str, isoforms: str) -> str:
    """
    Converts a .bed file to .gtf

    @type bed: str
    @param bed: path to .bed file
    @type isoforms: str
    @param isoforms: path to the isoforms file
    """

    global log

    gtf = f"{bed.split('.')[0]}.gtf"
    cmd = f"{Constants.ToolNames.BED2GTF} {bed} {isoforms} {gtf}"
    sh = shell(cmd)

    log.record(
        f"using {Constants.ToolNames.BED2GTF} from {Constants.Metadata.BED2GTF_METADATA} to convert bed to gtf"
    )
    log.record(sh)
    log.record(f"gff file written to {gtf}")

    return gtf


def bed_to_gff(bed: str, isoforms: str) -> str:
    """
    Converts a .bed file to .gff

    @type bed: str
    @param bed: path to .bed file
    @type isoforms: str
    @param isoforms: path to the isoforms file
    """

    global log

    gff = f"{bed.split('.')[0]}.gff"
    cmd = f"{Constants.ToolNames.BED2GFF} {bed} {isoforms} {gff}"
    sh = shell(cmd)

    log.record(
        f"using {Constants.ToolNames.BED2GFF} from {Constants.Metadata.BED2GFF_METADATA} to convert bed to gff"
    )
    log.record(sh)
    log.record(f"gff file written to {gff}")

    return gff


def run(args: argparse.Namespace) -> None:
    """
    The postoga runner function

    @type args: subprocess.Namespace
    @param args: defined arguments
    """

    global log

    log = Log(args.path, Constants.FileNames.LOG)
    log.start()
    log.record(f"postoga started!")
    log.record(f"running with arguments: {vars(args)}")

    # Get the toga table and write isoforms
    table = toga_table(args)
    isoforms = write_isoforms(args, table)

    # Filtering step (if any)
    if any([args.by_class, args.by_rel, args.threshold]):
        bed = filter_bed(args, table)
    else:
        bed = os.path.join(args.path, Constants.FileNames.BED)

    # Conversion step
    if args.to == "gtf":
        bed_to_gtf(bed, isoforms)
    elif args.to == "gff":
        bed_to_gff(bed, isoforms)

    log.close()

    return


def parser():
    """Argument parser for postoga"""
    app = argparse.ArgumentParser()
    app.add_argument(
        "-p", "--path", help="Path to TOGA results directory", required=True, type=str
    )
    app.add_argument(
        "-bc",
        "--by-class",
        help="Filter parameter to only include certain orthology classes (I, PI, UL, M, PM, L, UL)",
        required=False,
        type=str,
    )
    app.add_argument(
        "-brel",
        "--by-rel",
        help="Filter parameter to only include certain orthology relationships (o2o, o2m, m2m, m2m, o2z)",
        required=False,
        type=str,
    )
    app.add_argument(
        "-thold",
        "--threshold",
        help="Filter parameter to preserve orthology scores greater or equal a given threshold (0.0 - 1.0)",
        required=False,
        type=str,
    )
    app.add_argument(
        "-to",
        "--to",
        help="Specify the conversion format for .bed (query_annotation/filtered) file (gtf, gff3)",
        required=True,
        type=str,
        choices=["gtf", "gff"],
    )

    if len(sys.argv) < 2:
        app.print_help()
        sys.exit(0)
    args = app.parse_args()

    return args


def main():
    args = parser()
    run(args)


if __name__ == "__main__":
    main()
