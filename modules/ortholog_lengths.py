""" Calculate the distribution of ortholog lengths in TOGA results."""

from modules.utils import shell
import pandas as pd
from constants import Constants
from logger import Log

__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "jose.gonzalesdezavala1@unmsm.edu.pe"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.7.0-devel"


def noel_lengths(path: str, model: str) -> str:
    """
    Calculate the distribution of ortholog lengths in the resulting gtf/gff file using NOEL.

    @type path: str
    @param path: path to the results directory
    @type model: str
    @param model: the gtf/gff gene model
    """
    log = Log.connect(path, Constants.FileNames.LOG)

    cmd = f"{Constants.ToolNames.NOEL} -g {model} -o {path}/{Constants.FileNames.LENGTHS}"
    sh = shell(cmd)

    info = [
        f"using {Constants.ToolNames.NOEL} from {Constants.Metadata.NOEL_METADATA} to calculate ortholog lengths",
        sh,
        f"lengths file written to {path}/{Constants.FileNames.LENGTHS}",
    ]

    [log.record(i) for i in info]

    return Constants.FileNames.LENGTHS


def process_lenghts(path: str, lengths: str) -> pd.Series:
    """
    Processes the lengths file to get the ortholog lengths.

    @type path: str
    @param path: path to the results directory
    @type lengths: str
    @param lengths: path to the lengths file
    """
    log = Log.connect(path, Constants.FileNames.LOG)

    lengths = f"{path}/{lengths}"
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


def calculate_lengths(path:str, model:str) -> pd.Series:
    """
    Calculate the distribution of ortholog lengths in the resulting gtf/gff file.

    @type path: str
    @param path: path to the results directory
    @type lengths: str
    @param lengths: path to the lengths file
    """
    lengths = noel_lengths(path, model)
    result = process_lenghts(path, lengths)

    return result
