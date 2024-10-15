#!/usr/bin/env python3


"""A module to filter the original .bed file based on the query table."""

import pandas as pd
import polars as pl
import os
from constants import Constants
from logger import Log
from modules.utils import bed_reader
from typing import Tuple, List, Dict


__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "jose.gonzalesdezavala1@unmsm.edu.pe"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.9.3-devel"


def filter_bed(
    togadir: str | os.PathLike,
    outdir: str | os.PathLike,
    table: pd.DataFrame,
    by_class: str | None,
    by_rel: str | None,
    threshold: float | None,
    paralog: float | None,
    engine: str = "pandas",
) -> Tuple:
    """
    Filters the original .bed file to produce a custom filtered file

    Parameters
    ----------
        togadir : str | os.PathLike
            The path to the TOGA results directory.
        outdir : str | os.PathLike
            The path to the output directory.
        table : pd.DataFrame
            The query table.
        by_class : list
            The classes to filter the table by.
        by_rel : list
            The relationships to filter the table by.
        threshold : float
            The orthology score threshold.
        paralog : float
            The paralogy score threshold.

    Returns
    -------
        tuple
            The filtered table and the custom table.

    Example
    -------
        >>> from modules.filter_query_annotation import filter_bed
    """

    log = Log.connect(outdir, Constants.FileNames.LOG)
    f = os.path.join(outdir, Constants.FileNames.FILTERED_BED)
    initial = len(table)

    if threshold:
        if engine != "polars":
            table = table[table["pred"] >= float(threshold)]
        else:
            table = table.filter(pl.col("pred") >= threshold)
        log.record(
            f"discarded {initial - len(table)} projections with orthology scores <{threshold}"
        )

    if by_class:
        edge = len(table)
        if engine != "polars":
            table = table[table["class"].isin(by_class.split(","))]
        else:
            table = table.filter(pl.col("class").is_in(by_class.split(",")))
        log.record(
            f"discarded {edge - len(table)} projections with classes other than {by_class}"
        )

    if by_rel:
        edge = len(table)
        if engine != "polars":
            table = table[table["relation"].isin(by_rel.split(","))]
        else:
            table = table.filter(pl.col("relation").is_in(by_rel.split(",")))
        log.record(
            f"discarded {edge - len(table)} projections with relationships other than {by_rel}"
        )
    if paralog:
        edge = len(table)
        if engine != "polars":
            table = table.groupby("helper").filter(
                lambda x: (x["pred"] > float(paralog)).sum() <= 1
            )
        else:
            table = table.filter(
                pl.col("helper").is_in(
                    table.group_by("helper")
                    .agg(
                        [(pl.col("pred") > paralog).sum().alias("count_above_paralog")]
                    )
                    .filter(pl.col("count_above_paralog") <= 1)["helper"]
                )
            )

        log.record(
            f"discarded {edge - len(table)} transcripts with more than 1 chain with orthology probs >{paralog}"
        )

    # read the original .bed file and filter it based on the transcripts table
    if engine != "polars":
        bed = pd.read_csv(
            os.path.join(togadir, Constants.FileNames.BED), sep="\t", header=None
        )
        bed = bed[bed[3].isin(table["transcripts"])]
        custom_table = table[table["transcripts"].isin(bed[3])]

        bed.to_csv(
            f,
            sep="\t",
            header=None,
            index=False,
        )

        stats = [
            custom_table["class"].value_counts().to_dict(),
            custom_table["relation"].value_counts().to_dict(),
        ]
    else:
        bed = pl.read_csv(
            os.path.join(togadir, Constants.FileNames.BED),
            separator="\t",
            has_header=False,
        )
        bed = bed.filter(pl.col("column_4").is_in(table["transcripts"]))
        custom_table = table.filter(pl.col("transcripts").is_in(bed["column_4"]))

        bed.write_csv(f, include_header=False, separator="\t")

        stats = [
            dict(
                zip(
                    *table.group_by("class")
                    .agg(pl.len())
                    .to_dict(as_series=False)
                    .values()
                )
            ),
            dict(
                zip(
                    *table.group_by("relation")
                    .agg(pl.len())
                    .to_dict(as_series=False)
                    .values()
                )
            ),
        ]

    info = [
        f"kept {len(bed)} projections after filters, discarded {initial - len(bed)}.",
        f"{len(bed)} projections are coming from {len(custom_table['helper'].unique())} unique transcripts and {len(custom_table['t_gene'].unique())} genes",
        f"class stats of new bed: {custom_table['class'].value_counts().to_dict()}",
        f"relation stats of new bed: {custom_table['relation'].value_counts().to_dict()}",
        f"filtered bed file written to {f}",
    ]

    [log.record(i) for i in info]

    return f, stats, len(custom_table["t_gene"].unique()), custom_table


def get_stats_from_bed(
    bed_path: str, table: pd.DataFrame | pl.DataFrame, engine: str = "pandas"
) -> Tuple[List[Dict], int]:
    """
    Get the stats of a given bed file

    Parameters
    ----------
        bed_path : str
            The path to the bed file.
        table : pd.DataFrame | pl.DataFrame
            The table to get the stats from.
        engine : str
            The engine to use. Default is "pandas".

    Returns
    -------
        tuple
            The stats and the number of transcripts in the bed file.

    Example
    -------
        >>> from modules.filter_query_annotation import get_stats_from_bed
        >>> get_stats_from_bed("path/to/bed", table)
    """
    bed = bed_reader(bed_path)

    if engine != "polars":
        bed_table = table[table["transcripts"].isin(bed[3])]
        stats = [
            bed_table["class"].value_counts().to_dict(),
            bed_table["relation"].value_counts().to_dict(),
        ]
    else:
        bed_table = table.filter(pl.col("transcripts").is_in(bed["column_4"]))
        stats = [
            dict(
                zip(
                    *table.group_by("class")
                    .agg(pl.len())
                    .to_dict(as_series=False)
                    .values()
                )
            ),
            dict(
                zip(
                    *table.group_by("relation")
                    .agg(pl.len())
                    .to_dict(as_series=False)
                    .values()
                )
            ),
        ]

    return stats, len(bed_table["t_gene"].unique())
