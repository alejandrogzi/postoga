#!/usr/bin/env python3


"""A module to filter the original .bed file based on the query table."""

import os
from typing import Dict, List, Optional, Tuple, Union

import pandas as pd
import polars as pl

from constants import Constants
from logger import Log
from modules.utils import bed_reader

__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "jose.gonzalesdezavala1@unmsm.edu.pe"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.9.3-devel"


def filter_bed(
    togadir: Union[str, os.PathLike],
    outdir: Union[str, os.PathLike],
    table: Union[pd.DataFrame, pl.DataFrame],
    by_class: Optional[Union[str, os.PathLike]],
    by_rel: Optional[Union[str, os.PathLike]],
    threshold: Optional[float],
    paralog: Optional[float],
    engine: Union[str, os.PathLike] = "pandas",
) -> Tuple:
    """
    Filters the original .bed file to produce a custom filtered file

    Parameters
    ----------
        togadir : Union[str, os.PathLike]|Union[str, os.PathLike]
            The path to the TOGA results directory.
        outdir : Union[str, os.PathLike]|Union[str, os.PathLike]
            The path to the output directory.
        table : pd.DataFrame
            The query table.
        by_class : list
            The classes to filter the table by.
        by_rel : list
            The orthology_relationships to filter the table by.
        threshold : float
            The orthology orthology_probability threshold.
        paralog : float
            The paralogy orthology_probability threshold.

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
            table = table[table["orthology_probability"] >= float(threshold)]
        else:
            table = table.filter(pl.col("orthology_probability") >= threshold)
        log.record(
            f"discarded {initial - len(table)} projections with orthology orthology_probabilitys <{threshold}"
        )

    if by_class:
        edge = len(table)
        if engine != "polars":
            table = table[table["status"].isin(by_class.split(","))]
        else:
            table = table.filter(pl.col("status").is_in(by_class.split(",")))
        log.record(
            f"discarded {edge - len(table)} projections with classes other than {by_class}"
        )

    if by_rel:
        edge = len(table)
        if engine != "polars":
            table = table[table["orthology_relation"].isin(by_rel.split(","))]
        else:
            table = table.filter(pl.col("orthology_relation").is_in(by_rel.split(",")))
        log.record(
            f"discarded {edge - len(table)} projections with orthology_relationships other than {by_rel}"
        )
    if paralog:
        edge = len(table)
        if engine != "polars":
            table = table.groupby("reference_transcript").filter(
                lambda x: (x["orthology_probability"] > float(paralog)).sum() <= 1
            )
        else:
            table = table.filter(
                pl.col("reference_transcript").is_in(
                    table.group_by("reference_transcript")
                    .agg(
                        [
                            (pl.col("orthology_probability") > paralog)
                            .sum()
                            .alias("count_above_paralog")
                        ]
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
        bed = bed[bed[3].isin(table["projection"])]
        custom_table = table[table["projection"].isin(bed[3])]

        bed.to_csv(
            f,
            sep="\t",
            header=None,
            index=False,
        )

        stats = [
            custom_table["status"].value_counts().to_dict(),
            custom_table["orthology_relation"].value_counts().to_dict(),
        ]
    else:
        bed = pl.read_csv(
            os.path.join(togadir, Constants.FileNames.BED),
            separator="\t",
            has_header=False,
        )
        bed = bed.filter(pl.col("column_4").is_in(table["projection"]))
        custom_table = table.filter(pl.col("projection").is_in(bed["projection"]))

        bed.write_csv(f, include_header=False, separator="\t")

        stats = [
            dict(
                zip(
                    *table.group_by("status")
                    .agg(pl.len())
                    .to_dict(as_series=False)
                    .values()
                )
            ),
            dict(
                zip(
                    *table.group_by("orthology_relation")
                    .agg(pl.len())
                    .to_dict(as_series=False)
                    .values()
                )
            ),
        ]

    info = [
        f"kept {len(bed)} projections after filters, discarded {initial - len(bed)}.",
        f"{len(bed)} projections are coming from {len(custom_table['reference_transcript'].unique())} unique transcripts and {len(custom_table['reference_gene'].unique())} genes",
        f"class stats of new bed: {custom_table['status'].value_counts().to_dict()}",
        f"orthology_relation stats of new bed: {custom_table['orthology_relation'].value_counts().to_dict()}",
        f"filtered bed file written to {f}",
    ]

    [log.record(i) for i in info]

    return f, stats, len(custom_table["reference_gene"].unique()), custom_table


def get_stats_from_bed(
    bed_path: str,
    table: Union[pd.DataFrame, pl.DataFrame],
    engine: Union[str, os.PathLike] = "pandas",
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
        bed_table = table[table["projection"].isin(bed[3])]
        stats = [
            bed_table["status"].value_counts().to_dict(),
            bed_table["orthology_relation"].value_counts().to_dict(),
        ]
    else:
        bed_table = table.filter(pl.col("projection").is_in(bed["column_4"]))
        stats = [
            dict(
                zip(
                    *table.group_by("status")
                    .agg(pl.len())
                    .to_dict(as_series=False)
                    .values()
                )
            ),
            dict(
                zip(
                    *table.group_by("orthology_relation")
                    .agg(pl.len())
                    .to_dict(as_series=False)
                    .values()
                )
            ),
        ]

    return stats, len(bed_table["reference_gene"].unique())
