#!/usr/bin/env python3

"""A module with postoga base utility functions."""

import os
import subprocess
from typing import List, Tuple, Union

import pandas as pd
import polars as pl

from constants import Constants

pd.options.mode.chained_assignment = None  # default='warn'

__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "jose.gonzalesdezavala1@unmsm.edu.pe"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.9.3-devel"


def shell(cmd: str) -> str:
    """
    Run a shell command and return the output as a string

    Parameters
    ----------
        cmd : str
            The shell command to run.

    Returns
    -------
        str
            The output of the shell command.

    Example
    -------
        >>> from modules.utils import shell
        >>> shell("ls -l")
    """
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True, check=True)
    return result.stdout.strip()


def bed_reader(
    bed: Union[str, os.PathLike], engine: str = "pandas"
) -> Union[pd.DataFrame, pl.DataFrame]:
    """
    Reads a .bed file and returns a pandas DataFrame

    Parameters
    ----------
        bed : str | os.PathLike
            The path to the .bed file.
        engine : str
            The engine to use to read the .bed file.

    Returns
    -------
        pd.DataFrame
            The .bed file as a pandas DataFrame.

    Example
    -------
        >>> from modules.utils import bed_reader
        >>> bed_reader("path/to/file.bed", engine="pandas")
    """
    if engine != "polars":
        return pd.read_csv(bed, sep="\t", header=None)
    else:
        return pl.read_csv(bed, separator="\t", has_header=False)


def ancestral_reader(ancestral: str, source: str, engine: str = "pandas") -> List:
    """
    Reads an ancestral file and returns a pandas DataFrame

    Parameters
    ----------
        ancestral : str | os.PathLike
            The path to the ancestral file.
        source : str
            The source column to extract from the ancestral file.

    Returns
    -------
        list
            The source column as a list.

    Example
    -------
        >>> from modules.utils import ancestral_reader
        >>> ancestral_reader("path/to/file.bed", "source")
    """
    if engine != "polars":
        df = pd.read_csv(ancestral, sep="\t").loc[:, source].to_list()
    else:
        df = (
            pl.read_csv(ancestral, separator="\t")
            .select(source)
            .to_dict(as_series=False)
            .get(source)
        )

    return df


def isoform_writer(
    outdir: Union[str, os.PathLike],
    table: Union[pd.DataFrame, pl.DataFrame],
    engine: str = "pandas",
) -> Tuple[str, str]:
    """
    Writes all isoforms to a text file

    Parameters
    ----------
        outdir : str | os.PathLike
            The output directory.
        table : pd.DataFrame | pl.DataFrame
            The DataFrame with the isoforms.
        engine : str
            The engine to use to write the .txt file.

    Returns
    -------
        str
            The path to the .txt file.

    Example
    -------
        >>> from modules.write_isoforms import isoform_writer
        >>> isoform_writer("path/to/output/directory", table, engine="pandas")
    """

    f = os.path.join(outdir, Constants.FileNames.OWNED_ISOFORMS)

    # get only gene:transcript pairs
    if engine != "polars":
        table.iloc[:, [1, 3]].dropna().to_csv(f, sep="\t", header=None, index=False)
    else:
        table.select("t_gene", "projection").drop_nulls().write_csv(
            f, separator="\t", include_header=False, quote_char=""
        )

    return f, f"gene-to-projection hash with {len(table)} entries written to {f}"
