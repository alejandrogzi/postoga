#!/usr/bin/env python3

"""A module to handle the haplotype branch of postoga."""

import os
import pandas as pd
import numpy as np
from modules.utils import bed_reader
from modules.make_query_table import query_table
from constants import Constants
from functools import reduce
from logger import Log

__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "jose.gonzalesdezavala1@unmsm.edu.pe"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.9.3-devel"


def get_haplotype_classes(paths: list, source: str) -> None:
    """
    @type paths: str
    @param paths: comma separated paths to TOGA results directories
    @type source: str
    @param source: the source of the haplotype classes

    cmd: postoga haplotypes -hpath path1,path2,path3 --source [query, loss]
    """

    log = Log.connect(paths[0], Constants.FileNames.LOG)
    dfs = []

    if source != "loss":
        # For each path, build a query table and filter it based on the bed file, then append to dfs
        for path in paths:
            table = query_table(path)
            bed = bed_reader(os.path.join(path, Constants.FileNames.BED))
            df = table[table["transcripts"].isin(bed[3])]
            dfs.append(df)
    else:
        # For each path, read loss_summ_data.tsv and append to dfs
        for path in paths:
            df = pd.read_csv(
                os.path.join(path, Constants.FileNames.CLASS),
                sep="\t",
                header=None,
                names=["projection", "transcripts", "class"],
            )
            dfs.append(df)
            log.record(
                f"number of projections in {path}: {len(df[df['projection'] == 'PROJECTION'])}, transcripts: {len(df[df['projection'] == 'TRANSCRIPT'])}, genes: {len(df[df['projection'] == 'GENE'])}, total: {len(df)}"
            )

    return dfs


def get_haplotype_rules(rule: list) -> dict:
    """
    @type rule: str
    @param rule: a string with the rule to to keep the "best" class for a given gene
    @rtype: dict
    @return: a dictionary of the form {rule:order}

    cmd: postoga -m haplotypes -hpath path1,path2 -rule "I>PI>UL>L>M>PM>PG>abs"
    """

    # Build a dictionary with the rules and adds "abs" at the end
    rules = {char: i for i, char in enumerate(rule)}
    rules["NF"] = len(rules)

    return rules


def merge_haplotypes(paths: list, source: str, rule: list) -> pd.DataFrame:
    """
    @type paths: str
    @param paths: comma separated paths to TOGA results directories
    @type source: str
    @param source: the source of the haplotype classes
    @type rule: list
    @param rule: a list with the rule to keep the "best" class for a given gene

    cmd: postoga haplotypes -hpath path1,path2,path3 --source [query, loss] --rule "I>PI>UL>L>M>PM>PG>abs"
    """

    log = Log.connect(paths[0], Constants.FileNames.LOG)

    dfs = get_haplotype_classes(paths, source)
    rules = get_haplotype_rules(rule)

    if len(dfs) < 2:
        raise ValueError("You must provide at least two paths to merge haplotypes")
    elif len(dfs) > 2:
        count = 0

        if source != "loss":
            for df in dfs:
                df.rename(
                    columns={
                        "t_gene": f"t_gene_{count}",
                        "class": f"class_{count}",
                        "helper": f"helper_{count}",
                        "relation": f"relation_{count}",
                    },
                    inplace=True,
                )
                count += 1
        else:
            for df in dfs:
                df.rename(
                    columns={
                        "projection": f"projection_{count}",
                        "class": f"class_{count}",
                    },
                    inplace=True,
                )
                count += 1

        multiple_table = reduce(
            lambda x, y: pd.merge(x, y, on="transcripts", how="outer"), dfs
        )
        multiple_table["consensus"] = None
        multiple_table.fillna("NF", inplace=True)

        for index, row in multiple_table.iterrows():
            classes = [row[f"class_{i}"] for i in range(len(dfs))]

            # Calculate the consensus class based on the provided rules
            consensus_class = min(classes, key=lambda c: rules[c])
            multiple_table.at[index, "consensus"] = consensus_class

        if source != "loss":
            for i in range(len(dfs)):
                multiple_table[f"t_gene_{i}"].replace("NF", np.NaN, inplace=True)
                multiple_table[f"helper_{i}"].replace("NF", np.NaN, inplace=True)
                multiple_table[f"relation_{i}"].replace("NF", np.NaN, inplace=True)
                if i > 0:
                    multiple_table["t_gene_0"].fillna(
                        multiple_table[f"t_gene_{i}"], inplace=True
                    )
                    multiple_table["helper_0"].fillna(
                        multiple_table[f"helper_{i}"], inplace=True
                    )
                    multiple_table["relation_0"].fillna(
                        multiple_table[f"relation_{i}"], inplace=True
                    )

            multiple_table = multiple_table[
                [
                    "t_gene_0",
                    "helper_0",
                    "transcripts",
                    "relation_0",
                    "consensus",
                ]
            ]

            log.record(
                f"number of projections in merged table: {len(multiple_table['transcripts'])}, transcripts: {len(multiple_table['helper_0'].unique())}, genes: {len(multiple_table['t_gene_0'].unique())}, total: {len(multiple_table)}"
            )
            log.record(
                f"class stats in merged table: {multiple_table['consensus'].value_counts().to_dict()}"
            )

        else:
            for i in range(len(dfs)):
                multiple_table[f"projection_{i}"].replace("NF", np.NaN, inplace=True)
                if i > 0:
                    multiple_table["projection_0"].fillna(
                        multiple_table[f"projection_{i}"], inplace=True
                    )

            multiple_table = multiple_table[
                ["projection_0", "transcripts", "consensus"]
            ]

            log.record(
                f"number of projections in merged table: {len(multiple_table[multiple_table['projection_0'] == 'PROJECTION'])}, transcripts: {len(multiple_table[multiple_table['projection_0'] == 'TRANSCRIPT'])}, genes: {len(multiple_table[multiple_table['projection_0'] == 'GENE'])}, total: {len(multiple_table)}"
            )
            log.record(
                f"class stats in merged table: {multiple_table['consensus'].value_counts().to_dict()}"
            )

        multiple_table.to_csv(
            os.path.join(paths[0], Constants.FileNames.HAPLOTYPE),
            sep="\t",
            index=False,
            header=False,
        )

        log.record(
            f"consensus table written to { os.path.join(paths[0], Constants.FileNames.HAPLOTYPE)}"
        )

        return multiple_table

    else:
        paired_table = reduce(
            lambda x, y: pd.merge(x, y, on="transcripts", how="outer"), dfs
        )
        paired_table["consensus"] = None
        paired_table["class_x"].fillna("NF", inplace=True)
        paired_table["class_y"].fillna("NF", inplace=True)

        for index, row in paired_table.iterrows():
            x = row["class_x"]
            y = row["class_y"]

            if x == y:
                paired_table.at[index, "consensus"] = x
            elif rules[x] < rules[y]:
                paired_table.at[index, "consensus"] = x
            else:
                paired_table.at[index, "consensus"] = y

        if source != "loss":
            paired_table["t_gene_x"].fillna(paired_table["t_gene_y"], inplace=True)
            paired_table["helper_x"].fillna(paired_table["helper_y"], inplace=True)
            paired_table["relation_x"].fillna(paired_table["relation_y"], inplace=True)

            paired_table = paired_table[
                [
                    "t_gene_x",
                    "helper_x",
                    "transcripts",
                    "relation_x",
                    "consensus",
                ]
            ]

            log.record(
                f"number of projections in merged table: {len(paired_table['transcripts'])}, transcripts: {len(paired_table['helper_x'].unique())}, genes: {len(paired_table['t_gene_x'].unique())}, total: {len(paired_table)}"
            )
            log.record(
                f"class stats in merged table: {paired_table['consensus'].value_counts().to_dict()}"
            )

        else:
            paired_table["projection_x"].fillna(
                paired_table["projection_y"], inplace=True
            )
            paired_table = paired_table[["projection_x", "transcripts", "consensus"]]

            log.record(
                f"number of projections in merged table: {len(paired_table[paired_table['projection_x'] == 'PROJECTION'])}, transcripts: {len(paired_table[paired_table['projection_x'] == 'TRANSCRIPT'])}, genes: {len(paired_table[paired_table['projection_x'] == 'GENE'])}, total: {len(paired_table)}"
            )
            log.record(
                f"class stats in merged table: {paired_table['consensus'].value_counts().to_dict()}"
            )

        paired_table.to_csv(
            os.path.join(paths[0], Constants.FileNames.HAPLOTYPE), sep="\t", index=False
        )

        log.record(
            f"consensus table written to { os.path.join(paths[0], Constants.FileNames.HAPLOTYPE)}"
        )

        return paired_table
