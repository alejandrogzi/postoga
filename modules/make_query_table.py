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
from typing import Dict, Optional, Tuple, Union

import pandas as pd
import polars as pl

from constants import Constants
from logger import Log

__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "jose.gonzalesdezavala1@unmsm.edu.pe"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.9.3-devel"


MISSING_PLACEHOLDER = "GENE_NOT_FOUND"


def query_table(
    path: Union[str, os.PathLike],
    outdir: Union[str, os.PathLike],
    engine: str = "pandas",
) -> Union[pd.DataFrame, pl.DataFrame]:
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
    path: Union[str, os.PathLike],
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, Dict[str, str]]:
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
    score["tx_with_chain"] = score["transcript"] + "#" + score["chain"].astype(str)
    score = score[["tx_with_chain", "pred", "transcript"]]

    loss = pd.read_csv(os.path.join(path, Constants.FileNames.LOSS), sep="\t")

    inact_muts = pd.read_csv(
        os.path.join(path, Constants.FileNames.INACTIVATING_MUTATIONS), sep="\t"
    )

    query_genes = (
        pd.read_csv(os.path.join(path, Constants.FileNames.QUERY_GENES), sep="\t")
        .set_index("projection")
        .to_dict()["query_gene"]
    )

    return orthology, loss, score, inact_muts, query_genes


def make_pd_table(path: Union[str, os.PathLike]) -> pd.DataFrame:
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

    orthology, loss, score, inact_muts, query_genes = read_pd(path)

    # INFO: subsets loss to consider only projections
    loss = loss.query("level == 'PROJECTION'").copy()
    loss["helper"] = loss["entry"].str.rsplit(".", n=1).str[0]

    ortho_x_loss = pd.merge(
        orthology, loss, left_on="q_transcript", right_on="entry", how="outer"
    ).fillna({"helper": orthology["t_gene"]})

    table = pd.merge(
        ortho_x_loss, score, left_on="entry", right_on="tx_with_chain", how="outer"
    )
    table.fillna({"helper": table["transcript"]}, inplace=True)

    # INFO: create a new column with a rename orthology relationship
    table["relation"] = table["orthology_class"].map(Constants.ORTHOLOGY_TYPE)
    table.fillna({"tx_with_chain": table["transcript"]}, inplace=True)

    # INFO: creates a dictionary: transcript -> gene
    isoforms_dict = table.dropna().set_index("helper")["t_gene"].to_dict()
    table.fillna({"t_gene": table["helper"].map(isoforms_dict)}, inplace=True)

    # INFO: adds inactivating mutation data
    inact_muts = inact_muts.groupby("projection", as_index=False).agg(
        {"mut_id": lambda x: list(x)}
    )
    inact_muts["inact_mut"] = [True for x in range(len(inact_muts))]
    table = pd.merge(
        table, inact_muts, left_on="tx_with_chain", right_on="projection", how="outer"
    )

    table.fillna({"inact_mut": False, "pred": 0.0}, inplace=True)

    # INFO: filling query genes
    table.fillna({"q_gene": table["tx_with_chain"].map(query_genes)}, inplace=True)
    table.fillna({"q_gene": table["t_gene"]}, inplace=True)
    table.fillna({"q_gene": MISSING_PLACEHOLDER}, inplace=True)

    # INFO: completing projection names
    table.fillna({"tx_with_chain": table["entry"]}, inplace=True)

    info = [
        f"found {len(table)} projections, {len(table['helper'].unique())} unique transcripts, {len(table['t_gene'].unique())} unique genes",
        f"class stats: {table['status'].value_counts().to_dict()}",
        f"relation stats: {table['relation'].value_counts().to_dict()}",
        # f"joined fragments predictions: {len(table[table['chain'] == 1])}",
    ]

    table = table[
        [
            "t_gene",
            "q_gene",
            "helper",
            "tx_with_chain",
            "relation",
            "status",
            "pred",
            "inact_mut",
            "mut_id",
        ]
    ].rename(
        columns={
            "t_gene": "reference_gene",
            "q_gene": "query_gene",
            "helper": "reference_transcript",
            "tx_with_chain": "projection",
            "relation": "orthology_relation",
            "status": "status",
            "pred": "orthology_probability",
            "inact_mut": "has_inact_mut",
            "mut_id": "inact_mut_id",
        }
    )

    [log.record(i) for i in info]

    return table


def read_pl(
    path: Union[str, os.PathLike],
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
        (score["transcript"] + "#" + score["chain"]).alias("tx_with_chain")
    )["tx_with_chain", "pred", "transcript"]

    loss = pl.read_csv(
        os.path.join(path, Constants.FileNames.LOSS),
        separator="\t",
    )

    inact_muts = pl.read_csv(
        os.path.join(path, Constants.FileNames.INACTIVATING_MUTATIONS),
        separator="\t",
    )

    # WARN: dropping paralogs for now!
    # paras = pl.read_csv(
    #     os.path.join(path, Constants.FileNames.PARALOGS),
    #     separator="\t",
    #     has_header=False,
    #     new_columns=["transcripts"],
    # )
    # paralogs = (
    #     score.join(paras, on="transcripts")
    #     .group_by("gene")
    #     .agg(pl.col("pred").max())
    #     .rename({"gene": "helper", "pred": "paralog_prob"})
    # )

    return orthology, loss, score, inact_muts


def make_pl_table(path: Union[str, os.PathLike]) -> pl.DataFrame:
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
    orthology, loss, score, inact_muts = read_pl(path)

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
        "status",
        "pred",
        "q_gene",
        "paralog_prob",
        "chain",
    )
