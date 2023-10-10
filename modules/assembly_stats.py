#!/usr/bin/env python3

""" 
Re-implementation of TOGA_assemblyStats.py -stats script for postoga. 

assembly_stats.py is a module to benchmark assembly quality in
terms of completeness. This module is a re-implementation
of the original TOGA_assemblyStats.py script.

By default, this module calculates the assembly completeness
based on the Ancestral_placental.txt file.
"""


import pandas as pd
from constants import Constants
from logger import Log
from modules.utils import bed_reader, ancestral_reader


__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "jose.gonzalesdezavala1@unmsm.edu.pe"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.6.0-devel"


def get_classes(path: str, bed: str, table: pd.DataFrame) -> pd.DataFrame:
    """
    @type path: str
    @param path: path to the results directory
    @type bed: str
    @param bed: path to original/filtered bed file
    @type table: pd.DataFrame
    @param table: a pandas DataFrame
    """
    log = Log.connect(path, Constants.FileNames.LOG)

    # Creates a table with unique genes in the query annotation (base or filtered) and sort them based on their class
    bed = bed_reader(bed)
    genes = table[table["transcripts"].isin(bed[3])].sort_values(
        by=["t_gene", "class"], key=lambda x: x.map(Constants.ORDER)
    )

    genes.drop_duplicates("t_gene", keep="first", inplace=True)

    log.record(f"number of unique genes in query: {len(genes)}")

    return genes


def qual_by_ancestral(
    path: str, bed: str, table: pd.DataFrame, assembly_qual: str
) -> None:
    """
    @type path: str
    @param path: path to the results directory
    @type bed: str
    @param bed: path to original/filtered bed file
    @type table: pd.DataFrame
    @param table: a pandas DataFrame
    @type assembly_qual: str
    @param assembly_qual: path to the ancestral placental file
    """
    log = Log.connect(path, Constants.FileNames.LOG)

    genes = get_classes(path, bed, table)
    ancestral = ancestral_reader(assembly_qual)
    overlap = genes[genes["t_gene"].isin(ancestral)]

    stats = overlap["class"].value_counts().to_dict()

    log.record(
        f"number of ancestral genes/custom genes ({len(ancestral)}) in query: {len(overlap)}, {len(overlap)/len(ancestral)*100:.2f}% overlap"
    )
    log.record(f"ancestral genes in query by class: {stats}")

    return stats
