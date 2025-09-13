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
from typing import Dict, Tuple, Union

import pandas as pd
import numpy as np

from constants import Constants
from logger import Log

__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "jose.gonzalesdezavala1@unmsm.edu.pe"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.9.3-devel"


MISSING_PLACEHOLDER = "GENE_NOT_FOUND"


def query_table(
    path: Union[str, os.PathLike],
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
    return make_pd_table(path)


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


def make_pd_table(
    path: Union[str, os.PathLike],
) -> pd.DataFrame:
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
    table.fillna({"q_gene": table["entry"].map(query_genes)}, inplace=True)
    table.fillna({"q_gene": table["t_gene"]}, inplace=True)

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

    # INFO: patch on retro labels!
    table.reference_transcript = [
        tx.split("#")[0].split(".")[0] if tx is not np.nan else tx
        for tx in table.reference_transcript
    ]
    transcript_to_gene = (
        table.dropna(subset=["reference_gene"])
        .set_index("reference_transcript")["reference_gene"]
        .to_dict()
    )

    table["reference_gene"] = table["reference_gene"].fillna(
        table["reference_transcript"].map(transcript_to_gene)
    )
    table["query_gene"] = table["query_gene"].fillna(
        table["reference_transcript"].map(transcript_to_gene)
    )

    # INFO: reading query_annotation.bed!
    try:
        bed = pd.read_csv(
            os.path.join(path, Constants.FileNames.BED), sep="\t", header=None
        )
    except:
        bed = pd.read_csv(
            os.path.join(path, Constants.FileNames.BED_UTR), sep="\t", header=None
        )
    transition = bed[[3]].rename(columns={3: "projection"})

    transition["reference_transcript"] = [
        p.split("#")[0].split(".")[0] for p in transition.projection
    ]
    transition["query_gene"] = table["reference_transcript"].map(transcript_to_gene)
    transition["query_gene"] = [
        transition.iloc[idx, 2] + "#RETRO"
        if transition.iloc[idx, 0].split("#")[-1] == "retro"
        and transition.iloc[idx, 2] is not np.nan
        else transition.iloc[idx, 2]
        for idx in transition.index
    ]

    table = pd.concat(
        [table, transition[~transition.projection.isin(table.projection)]]
    )
    table.fillna({"query_gene": MISSING_PLACEHOLDER}, inplace=True)

    [log.record(i) for i in info]

    return table
