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
import polars as pl
from constants import Constants
from logger import Log
from modules.utils import bed_reader, ancestral_reader
from typing import Dict, List, Tuple


__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "jose.gonzalesdezavala1@unmsm.edu.pe"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.9.3-devel"


def qual_by_ancestral(
    outdir: str | os.PathLike,
    bed: str | os.PathLike,
    table: pd.DataFrame | pl.DataFrame,
    assembly_qual: str,
    source: str,
    engine: str = "pandas",
) -> Dict[str, int]:
    """
    Calculate the number of ancestral genes in the query annotation

    Parameters
    ----------
    outdir : str | os.PathLike
        The path to the output directory.
    bed : str | os.PathLike
        The path to the original/filtered bed file.
    table : pd.DataFrame | pl.DataFrame
        The query annotation as a pandas DataFrame.
    assembly_qual : str
        The path to the ancestral placental file.
    source : str
        The source of the ancestral placental file.
    engine : str
        The engine to use for data manipulation (default: "pandas").

    Returns
    -------
    None

    Example
    -------
    >>> from modules.assembly_stats import qual_by_ancestral
    >>> qual_by_ancestral("output", "data.bed", table, "ancestral.txt", "source")
    """
    log = Log.connect(outdir, Constants.FileNames.LOG)

    ancestral = ancestral_reader(assembly_qual, source, engine)
    bed = bed_reader(bed, engine)

    if engine != "polars":
        genes = (
            table[table["transcripts"].isin(bed[3])]
            .sort_values(by=["t_gene", "class"], key=lambda x: x.map(Constants.ORDER))
            .drop_duplicates("t_gene", keep="first")
        )

        overlap = genes[genes["t_gene"].isin(ancestral)]
        stats = overlap["class"].value_counts().to_dict()
    else:
        genes = table.filter(pl.col("transcripts").is_in(bed["column_4"])).unique(
            subset="t_gene", keep="first"
        )

        overlap = genes.filter(pl.col("t_gene").is_in(ancestral))
        stats = dict(
            zip(*table.group_by("class").len().to_dict(as_series=False).values())
        )

    log.record(
        f"number of ancestral genes/custom genes ({len(ancestral)}) in query: {len(overlap)}, {len(overlap)/len(ancestral)*100:.2f}% overlap"
    )
    log.record(f"ancestral genes in query by class: {stats}")

    return stats


def overlap_busco(
    outdir: str | os.PathLike,
    db: List[str],
    table: pd.DataFrame | pl.DataFrame,
    src: str,
    engine: str = "pandas",
) -> List[Tuple[str, float]]:
    """
    Calculate the overlap between the query annotation and curated BUSCO databases

    Parameters
    ----------
    outdir : str | os.PathLike
        The path to the output directory.
    db : List[str]
        The list of curated BUSCO databases.
    table : pd.DataFrame | pl.DataFrame
        The query annotation as a pandas DataFrame.
    src : str
        The source of the query annotation (ensembl, entrez, gene_name).
    engine : str
        The engine to use for data manipulation (default: "pandas").

    Returns
    -------
    List[Tuple[str, float]]
        Tuple with the database name and the percentage of overlap.

    Example
    -------
    >>> from modules.assembly_stats import overlap_busco
    >>> overlap_busco("output", ["mammals", "birds"], table, "ensembl")
    """

    log = Log.connect(outdir, Constants.FileNames.LOG)
    track = []

    for odb in db:
        log.record(f"running pseudo-BUSCO for {odb} database")
        odb_path = f"{Constants.FileNames.SUPPLY_FOLDER}/odbs/{odb}.txt"

        if engine != "polars":
            odb_df = pd.read_csv(odb_path, sep="\t")
            odb_df = odb_df.loc[odb_df[src].dropna().index]
        else:
            odb_df = pl.read_csv(odb_path, separator="\t").filter(
                pl.col(src).drop_nulls()
            )

        log.record(
            f"number of genes in {odb} database with {src} nomenclature: {len(odb_df)}"
        )

        overlap = set(odb_df[src]).intersection(set(table["t_gene"]))
        track.append((odb, len(overlap) / len(odb_df) * 100))

    return track


def busco_completeness(
    outdir: str | os.PathLike,
    table: pd.DataFrame | pl.DataFrame,
    src: str,
    phylo: str,
    engine: str = "pandas",
) -> List:
    """
    Calculate the completeness of the assembly based on curated BUSCO databases

    Parameters
    ----------
    outdir : str | os.PathLike
        The path to the output directory.
    table : pd.DataFrame | pl.DataFrame
        The query annotation as a pandas DataFrame.
    src : str
        The source of the query annotation (ensembl, entrez, gene_name).
    phylo : str
        The phylogenetic group to use for the analysis (mammals, birds).

    Returns
    -------
    List[]

    """

    log = Log.connect(outdir, Constants.FileNames.LOG)

    # base functionality: compare assembly against eukaryota and vertebrata BUSCO databases
    base = overlap_busco(outdir, list(Constants.BUSCO_DBS_BASE.values()), table, src)

    if phylo == "mammals":
        supp = overlap_busco(
            outdir, list(Constants.BUSCO_DBS_MAMMALS.values()), table, src
        )
    elif phylo == "aves":
        supp = overlap_busco(
            outdir, list(Constants.BUSCO_DBS_BIRDS.values()), table, src
        )
    else:
        return base

    stats = base + supp
    log.record(f"pseudo-BUSCO stats: {stats}")

    return stats
