#!/usr/bin/env python3

import os
import argparse
import sys
from constants import Constants
from logger import Log
from modules.make_query_table import query_table
from modules.utils import isoform_writer
from modules.filter_query_annotation import filter_bed, get_stats_from_bed
from modules.assembly_stats import qual_by_ancestral, busco_completeness
from modules.haplotype_branch import merge_haplotypes
from modules.plotter import postoga_plotter
from modules.ortholog_lengths import calculate_lengths
from typing import Type
from rustools import convert, extract_seqs

__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "jose.gonzalesdezavala1@unmsm.edu.pe"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.9.3-devel"


class TogaDir:
    """A class to represent a TOGA results directory."""

    def __init__(self, args: argparse.Namespace) -> None:
        """
        Constructs all the necessary attributes for the TogaDir object.

        Parameters:
        ----------
        args : argparse.Namespace
            The argparse namespace object with all the arguments.

        Returns:
        ----------
        None
        """

        self.args = args
        self.mode = args.mode
        self.outdir = os.path.abspath(args.outdir)

        if args.mode != "haplotype":
            """The default branch of postoga"""
            self.togadir = args.togadir
            self.to = args.to
            self.log = Log(self.outdir, Constants.FileNames.LOG)

            self.by_class = args.by_class if args.by_class else None
            self.by_rel = args.by_rel if args.by_rel else None
            self.threshold = args.threshold if args.threshold else None
            self.para_threshold = args.paralog if args.paralog else None

            self.engine = args.engine
            self.q_assembly = args.assembly_qual
            self.species = args.species
            self.source = args.source
            self.phylo = args.phylo
            self.plot = args.plot
            self.isoforms = args.isoforms
            self.extract = args.extract

            self.codon = os.path.join(self.togadir, Constants.FileNames.CODON)
            self.protein = os.path.join(self.togadir, Constants.FileNames.PROTEIN)
            self.filtered_codon = os.path.join(
                self.outdir, Constants.FileNames.FILTERED_CODON
            )
            self.filtered_protein = os.path.join(
                self.outdir, Constants.FileNames.FILTERED_PROTEIN
            )
        else:
            """The haplotype branch of postoga"""
            self.togadirs = args.haplotype_path.split(",")
            self.rule = args.rule.split(">")
            self.source = args.source
            self.log = Log(self.outdir, Constants.FileNames.LOG)

    def run(self) -> None:
        """
        Run the postoga pipeline.

        Returns:
        ----------
        None

        Examples:
        ----------
        >>> TogaDir(args).run()
        """

        os.mkdir(self.outdir) if not os.path.isdir(self.outdir) else sys.exit(
            f"Error: {self.outdir} already exists!"
        )

        self.log.start()
        self.log.intro()
        self.log.record(f"postoga started!")
        self.log.record(
            f"running in mode {self.mode} with arguments: {vars(self.args)}"
        )

        if self.mode != "haplotype":
            self.table = query_table(self.togadir, self.outdir, self.engine)

            if not self.isoforms:
                self.isoforms, msg = isoform_writer(self.outdir, self.table)
                self.log.record(msg)
            else:
                self.log.record(
                    f"using custom isoform table provided by the user: {self.isoforms}"
                )

            if any([self.by_class, self.by_rel, self.threshold, self.para_threshold]):
                self.bed, self.stats, self.ngenes, self.custom_table = filter_bed(
                    self.togadir,
                    self.outdir,
                    self.table,
                    self.by_class,
                    self.by_rel,
                    self.threshold,
                    self.para_threshold,
                    self.engine,
                )
                self.base_stats, _ = get_stats_from_bed(
                    os.path.join(self.togadir, Constants.FileNames.BED),
                    self.table,
                    self.engine,
                )
            else:
                self.bed = os.path.join(self.togadir, Constants.FileNames.BED)
                self.base_stats, self.ngenes = get_stats_from_bed(
                    self.bed, self.table, self.engine
                )
                self.stats = None
                self.custom_table = self.table

            if self.to == "gtf":
                self.gtf = os.path.join(
                    self.outdir,
                    f"{os.path.splitext(os.path.basename(self.bed))[0]}.gtf",
                )
                self.gmodel = convert(self.bed, self.gtf, self.isoforms)
                self.log.record(f"Coversion to GTF file completed! {self.gtf} created!")
            elif self.to == "gff":
                self.gff = os.path.join(
                    self.outdir,
                    f"{os.path.splitext(os.path.basename(self.bed))[0]}.gff",
                )
                self.gmodel = convert(self.bed, self.gff, self.isoforms)
                self.log.record(f"Coversion to GFF file completed! {self.gff} created!")
            elif self.to == "bed":
                self.gmodel = self.bed
                self.log.record(
                    f"Skipping conversion to GTF or GFF file. Filtering only the .bed file and writing the results to {self.outdir}!"
                )

            self.ancestral_stats = qual_by_ancestral(
                self.outdir,
                self.bed,
                self.custom_table,
                self.q_assembly,
                self.source,
                self.engine,
            )

            self.completeness_stats = busco_completeness(
                self.outdir, self.custom_table, self.source, self.phylo, self.engine
            )

            if self.extract:
                extract_seqs(self.bed, self.protein, output=self.filtered_protein)
                extract_seqs(self.bed, self.codon, output=self.filtered_codon)

                self.log.record(
                    f"Extracted QUERY sequences from your filtered .bed file to {self.filtered_protein} and {self.filtered_codon}"
                )
            elif self.extract == "reference":
                extract_seqs(
                    self.bed, self.protein, self.extract, self.filtered_protein
                )
                extract_seqs(self.bed, self.codon, self.extract, self.filtered_codon)

                self.log.record(
                    f"Extracted REFERENCE sequences from your filtered .bed file to {self.filtered_protein} and {self.filtered_codon}"
                )
            else:
                raise ValueError(
                    f"Error: {self.extract} is not a valid option for the --extract argument. Please choose either 'query' or 'reference'."
                )

            if not self.plot:
                self.log.record("skipping plotting and only filtering the .bed file")
            # else:
            #     postoga_plotter(
            #         self.outdir,
            #         self.table,
            #         self.ancestral_stats,
            #         self.base_stats,
            #         self.ngenes,
            #         self.ortholog_lengths,
            #         self.source,
            #         self.completeness_stats,
            #         self.species,
            #         self.stats,
            #     )

            self.log.close()

        else:
            hap_classes = merge_haplotypes(self.togadirs, self.source, self.rule)
            self.log.close()


def base_branch(subparsers, parent_parser):
    base_parser = subparsers.add_parser(
        "base", help="Base mode", parents=[parent_parser]
    )
    base_parser.add_argument(
        "--togadir",
        "-td",
        help="Path to TOGA results directory",
        required=True,
        type=str,
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
        default="gtf",
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
        "-p",
        "--plot",
        help="Flag to plot statistics about the filtered genes (default: False)",
        required=False,
        default=False,
        action="store_true",
    )
    base_parser.add_argument(
        "-par",
        "--paralog",
        help="Filter parameter to preserve transcripts with paralog projection probabilities less or equal to a given threshold (0.0 - 1.0)",
        required=False,
        type=float,
    )
    base_parser.add_argument(
        "-iso",
        "--isoforms",
        help="Path to a custom isoform table (default: None)",
        required=False,
        default=None,
        type=str,
    )
    base_parser.add_argument(
        "-e",
        "--engine",
        help="Database engine to create inner db representations (default: pandas)",
        required=False,
        choices=["pandas", "polars"],
        type=str,
        default="pandas",
    )
    base_parser.add_argument(
        "-ext",
        "--extract",
        help="Flag or option to extract sequences (only codon and protein alignments) from the filtered genes. "
        "Can be 'query', 'reference', or just set as a flag (default: False). When used as a flag extracting 'query' sequences is assumed.",
        required=False,
        default=False,
        nargs="?",
        const=True,
        choices=["query", "reference"],
    )


def haplotype_branch(subparsers, parent_parser):
    haplotype_parser = subparsers.add_parser(
        "haplotype", help="Haplotype mode", parents=[parent_parser]
    )
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


def parser() -> argparse.Namespace:
    """Argument parser for postoga"""
    app = argparse.ArgumentParser()
    parent_parser = argparse.ArgumentParser(add_help=False)

    parent_parser.add_argument(
        "--outdir",
        "-o",
        help="Path to posTOGA output directory",
        required=True,
        type=str,
    )

    subparsers = app.add_subparsers(dest="mode", help="Select mode")

    base_branch(subparsers, parent_parser)
    haplotype_branch(subparsers, parent_parser)

    if len(sys.argv) < 2:
        app.print_help()
        sys.exit(0)

    args = app.parse_args()

    return args
