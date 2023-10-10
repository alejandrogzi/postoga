#!/usr/bin/env python3


""" A moduule to write isoforms in gene:transcript pair format to a text file. """


import pandas as pd
import os
from constants import Constants
from logger import Log


__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "jose.gonzalesdezavala1@unmsm.edu.pe"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.6.0-devel"


def isoform_writer(path: str, table: pd.DataFrame) -> str:
    """
    Writes all isoforms to a text file

    @type path: str
    @param path: path to the results directory
    @type table: pd.DataFrame
    @param table: query table
    """

    log = Log.connect(path, Constants.FileNames.LOG)

    f = os.path.join(path, Constants.FileNames.OWNED_ISOFORMS)

    # Get only gene:transcript pairs
    table = table.iloc[:, [0, 2]]
    table.to_csv(f, sep="\t", header=None, index=False)

    log.record(f"gene-to-projection hash with {len(table)} entries written to {f}")

    return f
