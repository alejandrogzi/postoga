#!/usr/bin/env python3

""" A module with postoga base utility functions. """


import subprocess
import pandas as pd
import os
from constants import Constants

__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "jose.gonzalesdezavala1@unmsm.edu.pe"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.2.0-devel"


def shell(cmd: str) -> str:
    """
    Run a shell command and return the output as a string

    @type cmd: str
    @param cmd: shell command
    """
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return result.stdout.strip()


def bed_reader(bed: str) -> pd.DataFrame:
    """
    Reads a .bed file and returns a pandas DataFrame

    @type path: str
    @param path: path to .bed file
    """
    return pd.read_csv(bed, sep="\t", header=None)


def ancestral_reader(ancestral: str) -> list:
    """
    Reads an ancestral file and returns a pandas DataFrame

    @type ancestral: str
    @param ancestral: path to ancestral file
    """
    return pd.read_csv(ancestral, sep="\t", header=None).iloc[:, 0].to_list()
