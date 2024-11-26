""" Calculate the distribution of ortholog lengths in TOGA results."""

import os
from modules.utils import shell
import pandas as pd
from constants import Constants
from logger import Log
from typing import Union


__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "jose.gonzalesdezavala1@unmsm.edu.pe"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.7.0-devel"


def noel_lengths(outdir: Union[str, os.PathLike], model: str) -> str:
    """
    Calculate the distribution of ortholog lengths in the resulting gtf/gff file using NOEL.

    @type outdir: str | os.PathLike
    @param outdir: path to output directory
    @type model: str
    @param model: the gtf/gff gene model
    """
    log = Log.connect(outdir, Constants.FileNames.LOG)
    lengths = os.path.join(outdir, Constants.FileNames.LENGTHS)

    cmd = (
        f"{Constants.ToolNames.NOEL} -g {model} -o {lengths}"
    )
    sh = shell(cmd)

    info = [
        f"using {Constants.ToolNames.NOEL} from {Constants.Metadata.NOEL_METADATA} to calculate ortholog lengths",
        f"Running {cmd}",
        sh,
        f"lengths file written to {lengths}",
    ]

    [log.record(i) for i in info]

    return lengths


def process_lenghts(outdir: Union[str, os.PathLike], lengths: str) -> pd.Series:
    """
    Processes the lengths file to get the ortholog lengths.

    @type outdir: str | os.PathLike
    @param outdir: path to output directory
    @type lengths: str
    @param lengths: path to the lengths file
    """
    log = Log.connect(outdir, Constants.FileNames.LOG)

    df = pd.read_csv(lengths, sep="\t", header=None, usecols=[1], names=["lengths"])
    igenes = len(df)

    df = df[df["lengths"] < 10000]
    fgenes = len(df)
    ogenes = igenes - fgenes

    info = [
        f"processing lengths file",
        f"total genes: {igenes}",
        f"genes with lengths < 10000: {fgenes}",
        f"genes with lengths >= 10000: {ogenes}",
    ]

    [log.record(i) for i in info]

    return df["lengths"]


def calculate_lengths(outdir: Union[str, os.PathLike], model: str) -> pd.Series:
    """
    Calculate the distribution of ortholog lengths in the resulting gtf/gff file.

    @type outdir: str | os.PathLike
    @param outdir: path to output directory
    @type lengths: str
    @param lengths: path to the lengths file
    """
    lengths = noel_lengths(outdir, model)
    result = process_lenghts(outdir, lengths)

    return result
