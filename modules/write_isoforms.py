#!/usr/bin/env python3


""" A module to write isoforms in gene:transcript pair format to a text file. """


import pandas as pd
import os
from constants import Constants
from logger import Log

pd.options.mode.chained_assignment = None  # default='warn'

__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "jose.gonzalesdezavala1@unmsm.edu.pe"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.7.0-devel"


def isoform_writer(outdir: str | os.PathLike, table: pd.DataFrame) -> str:
    """
    Writes all isoforms to a text file

    @type outdir: str | os.PathLike
    @param outdir: path to output directory
    @type table: pd.DataFrame
    @param table: query table
    """

    log = Log.connect(outdir, Constants.FileNames.LOG)

    f = os.path.join(outdir, Constants.FileNames.OWNED_ISOFORMS)

    # Get only gene:transcript pairs
    table = table.iloc[:, [0, 2]]
    table.dropna(inplace=True)
    table.to_csv(f, sep="\t", header=None, index=False)

    log.record(f"gene-to-projection hash with {len(table)} entries written to {f}")

    return f
