#!/usr/bin/env python3


"""
Master script for postoga. 

postoga is a tool that automates the post-processing of TOGA results.
At its core, this tool takes a TOGA results directory and produces
a series of steps to reduce the amount of manual work required to pre-process
files for downstream analysis.
"""


import os
import argparse
import sys
from constants import Constants
from logger import Log
from modules.convert_from_bed import bed_to_gtf, bed_to_gff
from modules.make_query_table import query_table
from modules.write_isoforms import isoform_writer
from modules.filter_query_annotation import filter_bed
from modules.assembly_stats import qual_by_ancestral
from modules.haplotype_branch import merge_haplotypes


__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "jose.gonzalesdezavala1@unmsm.edu.pe"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.5.0-devel"


class TogaDir:
    """A class to represent a TOGA results directory."""

    def __init__(self, args: argparse.Namespace) -> None:
        """
        Constructs all the necessary attributes for the TogaDir object.

        @type args: argparse.Namespace
        @param args: defined arguments
        """

        self.mode = args.mode
        self.args = args

        if args.mode != "haplotype":
            """The default branch of postoga"""
            ##### STEP 1 #####
            self.path = args.path
            self.to = args.to
            self.q_assembly = args.assembly_qual
            self.log = Log(args.path, Constants.FileNames.LOG)
            self.by_class = args.by_class if args.by_class else None
            self.by_rel = args.by_rel if args.by_rel else None
            self.threshold = args.threshold if args.threshold else None
        else:
            """The haplotype branch of postoga"""
            self.paths = args.haplotype_path.split(",")
            self.rule = args.rule.split(">")
            self.source = args.source
            self.log = Log(self.paths[0], Constants.FileNames.LOG)

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
                self.bed = filter_bed(
                    self.path, self.table, self.by_class, self.by_rel, self.threshold
                )
            else:
                self.bed = os.path.join(self.path, Constants.FileNames.BED)

            if self.to == "gtf":
                self.gtf = bed_to_gtf(self.path, self.bed, self.isoforms)
            elif self.to == "gff":
                self.gff = bed_to_gff(self.path, self.bed, self.isoforms)

            ##### STEP 2 #####
            _ = qual_by_ancestral(self.path, self.bed, self.table, self.q_assembly)

            self.log.close()

        else:
            hap_classes = merge_haplotypes(self.paths, self.source, self.rule)
            self.log.close()


def base_branch(subparsers):
    base_parser = subparsers.add_parser("base", help="Base mode")
    base_parser.add_argument(
        "-p", "--path", help="Path to TOGA results directory", required=True, type=str
    )
    base_parser.add_argument(
        "-bc",
        "--by-class",
        help="Filter parameter to only include certain orthology classes (I, PI, UL, M, PM, L, UL)",
        required=False,
        type=str,
    )
    base_parser.add_argument(
        "-br",
        "--by-rel",
        help="Filter parameter to only include certain orthology relationships (o2o, o2m, m2m, m2m, o2z)",
        required=False,
        type=str,
    )
    base_parser.add_argument(
        "-th",
        "--threshold",
        help="Filter parameter to preserve orthology scores greater or equal to a given threshold (0.0 - 1.0)",
        required=False,
        type=str,
    )
    base_parser.add_argument(
        "-to",
        "--to",
        help="Specify the conversion format for .bed (query_annotation/filtered) file (gtf, gff3)",
        required=True,
        type=str,
        choices=["gtf", "gff"],
    )
    base_parser.add_argument(
        "-aq",
        "--assembly_qual",
        help="Calculate assembly quality based on a list of genes provided by the user (default: Ancestral_placental.txt)",
        required=False,
        type=str,
        default=Constants.FileNames.ANCESTRAL,
    )


def haplotype_branch(subparsers):
    haplotype_parser = subparsers.add_parser("haplotype", help="Haplotype mode")
    haplotype_parser.add_argument(
        "-hp",
        "--haplotype_path",
        help="Path to TOGA results directories separated by commas (path1,path2,path3)",
        required=True,
        type=str,
    )
    haplotype_parser.add_argument(
        "-r",
        "--rule",
        help="Rule to merge haplotype assemblies (default: I>PI>UL>L>M>PM>PG>abs)",
        required=False,
        type=str,
        default="I>PI>UL>L>M>PM>PG>NF",
    )
    haplotype_parser.add_argument(
        "-s",
        "--source",
        help="Source of the haplotype classes (query, loss)",
        required=False,
        type=str,
        choices=["query", "loss"],
        default="loss",
    )


def parser():
    """Argument parser for postoga"""
    app = argparse.ArgumentParser()
    subparsers = app.add_subparsers(dest="mode", help="Select mode")

    base_branch(subparsers)
    haplotype_branch(subparsers)

    if len(sys.argv) < 2:
        app.print_help()
        sys.exit(0)

    args = app.parse_args()

    return args


def main():
    args = parser()
    master = TogaDir(args)
    master.run()


if __name__ == "__main__":
    main()
