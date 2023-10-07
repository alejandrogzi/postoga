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
from modules.utils import shell, bed_reader
from modules.convert_from_bed import bed_to_gtf, bed_to_gff
from modules.make_query_table import query_table
from modules.write_isoforms import isoform_writer
from modules.assembly_stats import qual_by_ancestral

__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "jose.gonzalesdezavala1@unmsm.edu.pe"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.4.0-devel"


class TogaDir:
    """A class to represent a TOGA results directory."""

    def __init__(self, args: argparse.Namespace) -> None:
        """
        Constructs all the necessary attributes for the TogaDir object.

        @type args: argparse.Namespace
        @param args: defined arguments
        """

        self.log = Log(args.path, Constants.FileNames.LOG)
        self.mode = args.mode
        self.args = args

        if args.mode != "haplotype":
            """The default branch of postoga"""
            ##### STEP 1 #####
            self.path = args.path
            self.to = args.to
            self.q_assembly = args.assembly_qual

            if args.by_class:
                self.by_class = args.by_class
            if args.by_rel:
                self.by_rel = args.by_rel
            if args.threshold:
                self.threshold = args.threshold

        else:
            """The haplotype branch of postoga"""
            self.paths = args.haplotype_path.split(",")
            self.rule = args.rule.split(">")

    def filter_bed(self) -> str:
        """Filters the original .bed file to produce a custom filtered file"""

        initial = len(self.table)

        if self.threshold:
            table = self.table[self.table["pred"] >= float(self.threshold)]
            self.log.record(
                f"discarded {initial - len(table)} projections with orthology scores <{self.threshold}"
            )

        if self.by_class:
            edge = len(table)
            table = table[table["class"].isin(self.by_class.split(","))]
            self.log.record(
                f"discarded {edge - len(table)} projections with classes other than {self.by_class}"
            )

        if self.by_rel:
            edge = len(table)
            table = table[table["relation"].isin(self.by_rel.split(","))]
            self.log.record(
                f"discarded {edge - len(table)} projections with relationships other than {self.by_rel}"
            )

        # Read the original .bed file and filter it based on the transcripts table
        bed = pd.read_csv(
            os.path.join(self.path, Constants.FileNames.BED), sep="\t", header=None
        )
        bed = bed[bed[3].isin(table["transcripts"])]
        custom_table = table[table["transcripts"].isin(bed[3])]

        # Write the filtered .bed file
        f = os.path.join(self.path, Constants.FileNames.FILTERED_BED)
        bed.to_csv(
            f,
            sep="\t",
            header=None,
            index=False,
        )

        info = [
            f"kept {len(bed)} projections after filters, discarded {initial - len(bed)}.",
            f"{len(bed)} projections are coming from {len(custom_table['helper'].unique())} unique transcripts and {len(custom_table['t_gene'].unique())} genes",
            f"class stats of new bed: {custom_table['class'].value_counts().to_dict()}",
            f"relation stats of new bed: {custom_table['relation'].value_counts().to_dict()}",
            f"confidence stats of new bed: {custom_table['confidence_level'].value_counts().to_dict()}",
            f"filtered bed file written to {f}",
        ]

        [self.log.record(i) for i in info]

        return f

    def get_haplotype_classes(self) -> None:
        """
        @type paths: str
        @param paths: comma separated paths to TOGA results directories

        cmd: postoga -m haplotypes -hpath path1,path2,path3
        """
        paths = self.paths
        dfs = []

        # For each path build a query table and filtered based on queries annotated
        for path in paths:
            table = query_table(path)
            bed = bed_reader(path)
            df = table[table["transcripts"].isin(bed[3])]
            dfs.append(df)

        return dfs

    def run(self) -> None:
        """
        The postoga runner function

        @type args: subprocess.Namespace
        @param args: defined arguments
        """

        self.log.start()
        self.log.intro()
        self.log.record(f"postoga started!")
        self.log.record(
            f"running in mode {self.mode} with arguments: {vars(self.args)}"
        )

        if self.mode != "haplotype":
            self.table = query_table(self.path)
            self.isoforms = isoform_writer(self.path, self.table)

            if any([self.by_class, self.by_rel, self.threshold]):
                self.bed = self.filter_bed()
            else:
                self.bed = os.path.join(self.path, Constants.FileNames.BED)

            if self.to == "gtf":
                self.gtf = bed_to_gtf(self.path, self.bed, self.isoforms)
            elif self.to == "gff":
                self.gff = bed_to_gff(self.path, self.bed, self.isoforms)

            ##### STEP 2 #####
            qual_by_ancestral(self.path, self.bed, self.table, self.q_assembly)

        else:
            hap_classes = self.get_haplotype_classes()


def parser():
    """Argument parser for postoga"""
    app = argparse.ArgumentParser()
    app.add_argument(
        "-m",
        "--mode",
        help="Run mode",
        type=str,
        choices=["default", "haplotype"],
        default="default",
    )
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
    app.add_argument(
        "-aq",
        "--assembly_qual",
        help="Calculate assembly quality based on a list of genes provided by the user (default: Ancestral_placental.txt)",
        required=False,
        type=str,
        default=Constants.FileNames.ANCESTRAL,
    )
    app.add_argument(
        "-hpath",
        "--haplotype_path",
        help="Path to haplotype directory",
        required=False,
        type=str,
    )
    app.add_argument(
        "-r",
        "--rule",
        help="Rule to merge haplotype assemblies (default: I>PI>UL>L>M>PM>PG>abs)",
        required=False,
        type=str,
        default="I>PI>UL>L>M>PM>PG>abs",
    )

    if len(sys.argv) < 2:
        app.print_help()
        sys.exit(0)

    args = app.parse_args()

    if args.mode == "haplotype" and not args.haplotype_path:
        app.error(
            "haplotype mode requires the path to the haplotype directory, please provide it with the -hpath/--haplotype_path argument"
        )

    return args


def main():
    args = parser()
    master = TogaDir(args)
    master.run()


if __name__ == "__main__":
    main()
