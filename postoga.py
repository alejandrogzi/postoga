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
from modules.filter_query_annotation import filter_bed, get_stats_from_bed
from modules.assembly_stats import qual_by_ancestral, busco_completeness
from modules.haplotype_branch import merge_haplotypes
from modules.plotter import postoga_plotter
from modules.ortholog_lengths import calculate_lengths

__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "jose.gonzalesdezavala1@unmsm.edu.pe"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.8.0-devel"


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
        self.outdir = os.path.abspath(args.outdir)

        if args.mode != "haplotype":
            """The default branch of postoga"""
            ##### STEP 1 #####
            self.togadir = args.togadir
            self.to = args.to
            self.q_assembly = args.assembly_qual
            self.log = Log(self.outdir, Constants.FileNames.LOG)
            self.by_class = args.by_class if args.by_class else None
            self.by_rel = args.by_rel if args.by_rel else None
            self.threshold = args.threshold if args.threshold else None
            self.para_threshold = args.paralog if args.paralog else None
            self.species = args.species
            self.source = args.source
            self.phylo = args.phylo
            self.skip = args.skip
        else:
            """The haplotype branch of postoga"""
            self.togadirs = args.haplotype_path.split(",")
            self.rule = args.rule.split(">")
            self.source = args.source
            self.log = Log(self.outdir, Constants.FileNames.LOG)

    def run(self) -> None:
        """
        The postoga runner function

        @type args: subprocess.Namespace
        @param args: defined arguments
        """

        os.mkdir(self.outdir) if not os.path.isdir(self.outdir) else None
        self.log.start()
        self.log.intro()
        self.log.record(f"postoga started!")
        self.log.record(
            f"running in mode {self.mode} with arguments: {vars(self.args)}"
        )

        if self.mode != "haplotype":
            self.table = query_table(self.togadir)
            self.isoforms = isoform_writer(self.outdir, self.table)

            if any([self.by_class, self.by_rel, self.threshold, self.para_threshold]):
                self.bed, self.stats, self.ngenes, self.custom_table = filter_bed(
                    self.togadir, self.outdir, self.table, self.by_class, self.by_rel, self.threshold, self.para_threshold
                )
                self.base_stats, _ = get_stats_from_bed(
                    os.path.join(self.togadir, Constants.FileNames.BED), self.table
                )
            else:
                self.bed = os.path.join(self.togadir, Constants.FileNames.BED)
                self.base_stats, self.ngenes = get_stats_from_bed(self.bed, self.table)
                self.stats = None
                self.custom_table = self.table

            if self.to == "gtf":
                self.gmodel = bed_to_gtf(self.outdir, self.bed, self.isoforms)
            elif self.to == "gff":
                self.gmodel = bed_to_gff(self.outdir, self.bed, self.isoforms)
            elif self.to == "bed":
                self.gmodel = self.bed

            if self.skip:
                self.log.record("skipping steps 2, 3, and 4 and only filtering the .bed file")
            else:
                ##### STEP 2 #####
                self.ancestral_stats = qual_by_ancestral(
                    self.outdir, self.bed, self.custom_table, self.q_assembly, self.source
                )

                ##### STEP 3 #####
                self.ortholog_lengths = calculate_lengths(self.outdir, self.gmodel)

                ##### STEP 4 #####

                self.completeness_stats = busco_completeness(
                    self.outdir, self.custom_table, self.source, self.phylo
                )

                postoga_plotter(
                    self.togadir,
                    self.table,
                    self.ancestral_stats,
                    self.base_stats,
                    self.ngenes,
                    self.ortholog_lengths,
                    self.source,
                    self.completeness_stats,
                    self.species,
                    self.stats,
                )

            self.log.close()

        else:
            hap_classes = merge_haplotypes(self.togadirs, self.source, self.rule)
            self.log.close()


def base_branch(subparsers, parent_parser):
    base_parser = subparsers.add_parser("base", help="Base mode", parents = [parent_parser])
    base_parser.add_argument(
        "--togadir",
        "--td", 
        help="Path to TOGA results directory", required=True, type=str
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
        type=float,
    )
    base_parser.add_argument(
        "-to",
        "--to",
        help="Specify the conversion format for .bed (query_annotation/filtered) file (gtf, gff3) or just keep it as .bed (bed)",
        required=True,
        type=str,
        choices=["gtf", "gff", "bed"],
    )
    base_parser.add_argument(
        "-aq",
        "--assembly_qual",
        help="Calculate assembly quality based on a list of genes provided by the user (default: Ancestral_placental.txt)",
        required=False,
        type=str,
        default=Constants.FileNames.ANCESTRAL,
    )
    base_parser.add_argument(
        "-sp",
        "--species",
        help="Species name to be used as a reference for the assembly quality calculation (default: human)",
        required=False,
        choices=["human", "mouse", "chicken"],
        type=str,
        default=Constants.SPECIES_DEFAULT,
    )
    base_parser.add_argument(
        "-src",
        "--source",
        help="Source of the ancestral gene names (default: ENSG)",
        required=False,
        choices=["ensembl", "gene_name", "entrez"],
        type=str,
        default=Constants.SRC_DEFAULT,
    )
    base_parser.add_argument(
        "-phy",
        "--phylo",
        help="Phylogenetic group of your species (default: mammals)",
        required=False,
        choices=["mammals", "birds"],
        type=str,
        default=Constants.PHYLO_DEFAULT,
    )
    base_parser.add_argument(
        "-s",
        "--skip",
        help="Skip steps 2, 3, and 4 and only filter the .bed file",
        required=False,
        action="store_true",
    )
    base_parser.add_argument(
        "-par",
        "--paralog",
        help="Filter parameter to preserve transcripts with paralog projection probabilities less or equal to a given threshold (0.0 - 1.0)",
        required=False,
        type=float,
    )


def haplotype_branch(subparsers, parent_parser):
    haplotype_parser = subparsers.add_parser("haplotype", help="Haplotype mode", parents = [parent_parser])
    haplotype_parser.add_argument(
        "-hp",
        "--haplotype_dir",
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
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument(
        "--outdir",
        "-o",
        help="Path to posTOGA output directory",
        required=True,
        type=str
    )


    subparsers = app.add_subparsers(dest="mode", help="Select mode")


    base_branch(subparsers, parent_parser)
    haplotype_branch(subparsers, parent_parser)

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
