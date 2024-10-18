#!/usr/bin/env python3

"""
Module to create a query table directly from a TOGA output.

This module would expect the following the files within the
TOGA output directory:

.
├── meta
│   ├── classification_results
│   │   ├── orthology_scores.tsv
|   ├── paralogous_projections.tsv
│   └── orthology_resolution
│       └── orthology_classification.tsv
|
└── results
    ├── loss_summary.tsv
    └── orthology_classification.tsv

The inner constructor is written to /outdir/.toga.table
"""

import os
import pandas as pd
import polars as pl
from constants import Constants
from logger import Log
from typing import Tuple


__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "jose.gonzalesdezavala1@unmsm.edu.pe"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.9.3-devel"


def query_table(
    path: str | os.PathLike, outdir: str | os.PathLike, engine: str = "pandas"
) -> pd.DataFrame:
    """
    Constructs a query table from a TOGA output directory
    considering all possible projections and orthologs.

    Parameters:
    ----------
        path : str | os.PathLike

    Returns:
    ----------
        pd.DataFrame

    Examples:
    ----------
        >>> query_table("path/to/toga/output")
    """

    log = Log.connect(path, Constants.FileNames.LOG)

    if engine != "polars":
        table = make_pd_table(path)
        table.to_csv(
            os.path.join(outdir, Constants.FileNames.TOGA_TABLE), sep="\t", index=False
        )
    else:
        table = make_pl_table(path)
        table.write_csv(
            os.path.join(outdir, Constants.FileNames.TOGA_TABLE), separator="\t"
        )

    return table


def get_median(group: pd.DataFrame) -> float:
    """
    Get the median of a group

    Parameters:
    ----------
        group : pd.DataFrame

    Returns:
    ----------
        float

    Examples:
    ----------
        >>> get_median(group)
    """
    return group.loc[:, "pred"].median()


def read_pd(
    path: str | os.PathLike,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Pandas reader for TOGA directory

    Parameters:
    ----------
        path : str | os.PathLike

    Returns:
    ----------
        Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]

    Examples:
    ----------
        >>> read_pd("path/to/toga/output")
    """
    orthology = pd.read_csv(os.path.join(path, Constants.FileNames.ORTHOLOGY), sep="\t")

    score = pd.read_csv(os.path.join(path, Constants.FileNames.SCORES), sep="\t")
    score["transcripts"] = score["gene"] + "." + score["chain"].astype(str)
    score = score[["transcripts", "pred", "gene"]]

    loss = pd.read_csv(
        os.path.join(path, Constants.FileNames.LOSS),
        sep="\t",
        header=None,
        names=["projection", "transcript", "class"],
    )

    paras = pd.read_csv(
        os.path.join(path, Constants.FileNames.PARALOGS),
        sep="\t",
        header=None,
        names=["transcripts"],
    )
    paralogs = (
        pd.merge(score, paras, on="transcripts")
        .groupby("gene")
        .agg({"pred": "max"})
        .reset_index()
    )
    paralogs.columns = ["helper", "paralog_prob"]

    return orthology, loss, score, paralogs


def make_pd_table(path: str | os.PathLike) -> pd.DataFrame:
    """
    Constructs a query table from a TOGA output directory using pandas engine

    Parameters:
    ----------
        path : str | os.PathLike

    Returns:
    ----------
        pd.DataFrame

    Examples:
    ----------
        >>> make_pd_table("path/to/toga/output")
    """
    log = Log.connect(path, Constants.FileNames.LOG)

    orthology, loss, score, paralogs = read_pd(path)

    # subsets loss to consider only projections
    loss = loss[loss["projection"] == "PROJECTION"]
    loss["helper"] = loss["transcript"].str.rsplit(".", n=1).str[0]

    ortho_x_loss = pd.merge(
        orthology, loss, left_on="q_transcript", right_on="transcript", how="outer"
    ).fillna({"helper": orthology["t_gene"]})

    table = pd.merge(
        ortho_x_loss, score, left_on="transcript", right_on="transcripts", how="outer"
    )
    table.fillna({"helper": table["gene"]}, inplace=True)

    # create a new column with a rename orthology relationship
    table["relation"] = table["orthology_class"].map(Constants.ORTHOLOGY_TYPE)
    table.fillna({"transcripts": table["transcript"]}, inplace=True)

    # creates a dictionary: transcript -> gene
    isoforms_dict = table.dropna().set_index("helper")["t_gene"].to_dict()
    table.fillna({"t_gene": table["helper"].map(isoforms_dict)}, inplace=True)

    # add paralog probabilities
    table = pd.merge(table, paralogs, on="helper", how="left")
    table.fillna({"paralog_prob": 0}, inplace=True)

    # calculates median prediction values for joined fragments
    table["chain"] = [
        1 if x.split(".")[-1] == "-1" and y != "-" else 0
        for x, y in zip(table["transcripts"].fillna("-"), table["q_gene"].fillna("-"))
    ]
    medians = (
        table[(table["chain"] != 1) & (table["pred"] >= 0) & (table["q_gene"].isna())]
        .groupby("helper")
        .apply(get_median)
        .reset_index(name="npred")
    )
    table.loc[table["chain"] == 1, "pred"] = (
        table[table["chain"] == 1]
        .merge(medians, on="helper", how="left")["npred"]
        .to_list()
    )

    info = [
        f"found {len(table)} projections, {len(table['helper'].unique())} unique transcripts, {len(table['t_gene'].unique())} unique genes",
        f"class stats: {table['class'].value_counts().to_dict()}",
        f"relation stats: {table['relation'].value_counts().to_dict()}",
        f"joined fragments predictions: {len(table[table['chain'] == 1])}",
    ]

    [log.record(i) for i in info]

    return table[
        [
            "t_gene",
            "helper",
            "transcripts",
            "relation",
            "class",
            "pred",
            "q_gene",
            "paralog_prob",
            "chain",
        ]
    ]


def read_pl(
    path: str | os.PathLike,
) -> Tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    """
    Polars reader for TOGA directory

    Parameters:
    ----------
        path : str | os.PathLike

    Returns:
    ----------
        Tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame, pl.DataFrame]

    Examples:
    ----------
        >>> read_pl("path/to/toga/output")
    """
    orthology = pl.read_csv(
        os.path.join(path, Constants.FileNames.ORTHOLOGY), separator="\t"
    )

    score = pl.read_csv(os.path.join(path, Constants.FileNames.SCORES), separator="\t")
    score = score.with_columns(
        (score["gene"] + "." + score["chain"]).alias("transcripts")
    )["transcripts", "pred", "gene"]

    loss = pl.read_csv(
        os.path.join(path, Constants.FileNames.LOSS),
        separator="\t",
        has_header=False,
        new_columns=["projection", "transcript", "class"],
    )

    paras = pl.read_csv(
        os.path.join(path, Constants.FileNames.PARALOGS),
        separator="\t",
        has_header=False,
        new_columns=["transcripts"],
    )
    paralogs = (
        score.join(paras, on="transcripts")
        .group_by("gene")
        .agg(pl.col("pred").max())
        .rename({"gene": "helper", "pred": "paralog_prob"})
    )

    return orthology, loss, score, paralogs


def make_pl_table(path: str | os.PathLike) -> pl.DataFrame:
    """
    Constructs a query table from a TOGA output directory using polars engine

    Parameters:
    ----------
        path : str | os.PathLike

    Returns:
    ----------
        pl.DataFrame

    Examples:
    ----------
        >>> make_pl_table("path/to/toga/output")
    """
    log = Log.connect(path, Constants.FileNames.LOG)

    orthology, loss, score, paralogs = read_pl(path)

    loss = (
        loss.filter(pl.col("projection") == "PROJECTION")
        .with_columns(
            [
                pl.col("transcript").str.split(by=".").list.get(0).alias("helper1"),
                pl.col("transcript").str.split(by=".").list.get(1).alias("helper2"),
            ]
        )
        .with_columns((pl.col("helper1") + "." + pl.col("helper2")).alias("helper"))
        .drop("helper1", "helper2")
    )

    ortho_x_loss = orthology.join(
        loss, left_on="q_transcript", right_on="transcript", how="full"
    ).with_columns([pl.col("helper").fill_null(pl.col("t_gene"))])

    table = (
        ortho_x_loss.join(
            score, left_on="transcript", right_on="transcripts", how="full"
        )
        .with_columns([pl.col("helper").fill_null(pl.col("gene"))])
        .with_columns(
            [
                pl.col("orthology_class")
                .replace_strict(Constants.ORTHOLOGY_TYPE)
                .alias("relation"),
                pl.col("transcripts").fill_null(pl.col("transcript")),
            ]
        )
    )

    isoforms_dict = dict(
        zip(
            *table.select("helper", "t_gene")
            .drop_nulls()
            .unique()
            .to_dict(as_series=False)
            .values()
        )
    )

    table = (
        table.with_columns(
            pl.col("helper")
            .replace_strict(isoforms_dict, default="NOT FOUND")
            .alias("t_gene")
        )
        .join(paralogs, on="helper", how="left")
        .with_columns([pl.col("paralog_prob").fill_null(0)])
        .with_columns(
            pl.when(
                (
                    pl.col("transcripts").fill_null("-").str.split(".").arr.get(-1).str
                    == "-1"
                )
                & (pl.col("q_gene").fill_null("-") != "-")
            )
            .then(1)
            .otherwise(0)
            .alias("chain")
        )
    )

    medians = (
        table.filter(
            (pl.col("chain") != 1)
            & (pl.col("pred") >= 0)
            & (pl.col("q_gene").is_null())
        )
        .group_by("helper")
        .agg([pl.col("pred").median().alias("npred")])
    )

    # missing last step: integrating medians into the table when chain == 1

    return table.select(
        "t_gene",
        "helper",
        "transcripts",
        "relation",
        "class",
        "pred",
        "q_gene",
        "paralog_prob",
        "chain",
    )
