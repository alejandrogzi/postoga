#!/usr/bin/env python3

""" 
Re-implementation of TOGA_assemblyStats.py -stats script for postoga. 

assembly_stats.py is a module to benchmark assembly quality in
terms of completeness. This module is a re-implementation
of the original TOGA_assemblyStats.py script.

By default, this module calculates the assembly completeness
based on the Ancestral_placental.txt file.
"""

import os
import pandas as pd
from constants import Constants
from logger import Log
from modules.utils import bed_reader, ancestral_reader
from typing import Union


__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "jose.gonzalesdezavala1@unmsm.edu.pe"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.7.0-devel"


def get_classes(outdir: Union[str, os.PathLike], bed: str, table: pd.DataFrame) -> pd.DataFrame:
    """
    @type outdir: str | os.PathLike
    @param outdir: path to output directory
    @type bed: str
    @param bed: path to original/filtered bed file
    @type table: pd.DataFrame
    @param table: a pandas DataFrame
    """
    log = Log.connect(outdir, Constants.FileNames.LOG)

    # Creates a table with unique genes in the query annotation (base or filtered) and sort them based on their class
    bed = bed_reader(bed)
    genes = table[table["transcripts"].isin(bed[3])].sort_values(
        by=["t_gene", "class"], key=lambda x: x.map(Constants.ORDER)
    )

    genes.drop_duplicates("t_gene", keep="first", inplace=True)

    log.record(f"number of unique genes in query: {len(genes)}")

    return genes


def qual_by_ancestral(
    outdir: Union[str, os.PathLike], bed: str, table: pd.DataFrame, assembly_qual: str, source: str
) -> None:
    """
    @type outdir: str | os.PathLike
    @param outdir: path to output directory
    @type bed: str
    @param bed: path to original/filtered bed file
    @type table: pd.DataFrame
    @param table: a pandas DataFrame
    @type assembly_qual: str
    @param assembly_qual: path to the ancestral placental file
    """
    log = Log.connect(outdir, Constants.FileNames.LOG)

    genes = get_classes(outdir, bed, table)
    ancestral = ancestral_reader(assembly_qual, source)
    overlap = genes[genes["t_gene"].isin(ancestral)]

    stats = overlap["class"].value_counts().to_dict()

    log.record(
        f"number of ancestral genes/custom genes ({len(ancestral)}) in query: {len(overlap)}, {len(overlap)/len(ancestral)*100:.2f}% overlap"
    )
    log.record(f"ancestral genes in query by class: {stats}")

    return stats


def overlap_busco(outdir: Union[str, os.PathLike], db: str, table: pd.DataFrame, src: str) -> float:
    """
    @type outdir: str | os.PathLike
    @param outdir: path to output directory
    @type table: pd.DataFrame
    @param table: a pandas DataFrame
    @type src: str
    @param src: source of the query annotation (ensembl, entrez, gene_name)
    """

    log = Log.connect(outdir, Constants.FileNames.LOG)
    track = []

    for odb in db:
        log.record(f"running pseudo-BUSCO for {odb} database")
        odb_path = f"{Constants.FileNames.SUPPLY_FOLDER}/odbs/{odb}.txt"
        odb_df = pd.read_csv(odb_path, sep="\t")
        odb_df = odb_df.loc[odb_df[src].dropna().index]
        log.record(
            f"number of genes in {odb} database with {src} nomenclature: {len(odb_df)}"
        )

        overlap = set(odb_df[src]).intersection(set(table["t_gene"]))
        track.append((odb, len(overlap) / len(odb_df) * 100))

    return track


def busco_completeness(outdir: Union[str, os.PathLike], table: pd.DataFrame, src: str, phylo: str) -> list:
    """
    @type outdir: str | os.PathLike
    @param outdir: path to output directory
    @type table: pd.DataFrame
    @param table: a pandas DataFrame with the custom query table
    @type src: str
    @param src: source of the query annotation (ensembl, entrez, gene_name)
    @type phylo: str
    @param phylo: phylogenetic group of the query annotation (mammals, birds)
    """

    log = Log.connect(outdir, Constants.FileNames.LOG)

    # Base functionality: compare assembly against eukaryota and vertebrata BUSCO databases
    base = overlap_busco(outdir, Constants.BUSCO_DBS_BASE.values(), table, src)

    if phylo == "mammals":
        supp = overlap_busco(outdir, Constants.BUSCO_DBS_MAMMALS.values(), table, src)
    elif phylo == "aves":
        supp = overlap_busco(outdir, Constants.BUSCO_DBS_BIRDS.values(), table, src)

    stats = base + supp
    log.record(f"pseudo-BUSCO stats: {stats}")

    return stats
