#!/usr/bin/env python3


""" A module to convert .bed files to .gtf and .gff files. """


import os
from constants import Constants
from logger import Log
from modules.utils import shell


__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "jose.gonzalesdezavala1@unmsm.edu.pe"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.6.0-devel"


def bed_to_gtf(outdir: str | os.PathLike, bed: str, isoforms: str) -> str:
    """
    Converts a .bed file to .gtf

    @type outdir: str | os.PathLike,
    @param outdir: path to the output directory
    @type bed: str
    @param bed: path to .bed file
    @type isoforms: str
    @param isoforms: path to the isoforms file
    """

    log = Log.connect(outdir, Constants.FileNames.LOG)

    gtf = os.path.join(outdir, f"{os.path.splitext(os.path.basename(bed))[0]}.gtf")
    cmd = f"{Constants.ToolNames.BED2GTF} --bed {bed} --isoforms {isoforms} --output {gtf}"
    sh = shell(cmd)

    info = [
        f"using {Constants.ToolNames.BED2GTF} from {Constants.Metadata.BED2GTF_METADATA} to convert bed to gtf",
        sh,
        f"gtf file written to {gtf}",
    ]

    [log.record(i) for i in info]

    return gtf


def bed_to_gff(outdir: str | os.PathLike, bed: str, isoforms: str) -> str:
    """
    Converts a .bed file to .gff

    @type outdir: str | os.PathLike,
    @param outdir: path to the output directory
    @type bed: str
    @param bed: path to .bed file
    @type isoforms: str
    @param isoforms: path to the isoforms file
    """

    log = Log.connect(outdir, Constants.FileNames.LOG)

    gff = os.path.join(outdir, f"{os.path.splitext(os.path.basename(bed))[0]}.gff")
    cmd = f"{Constants.ToolNames.BED2GFF} --bed {bed} --isoforms {isoforms} --output {gff}"
    sh = shell(cmd)

    info = [
        f"using {Constants.ToolNames.BED2GFF} from {Constants.Metadata.BED2GFF_METADATA} to convert bed to gff",
        f"running {cmd}",
        sh,
        f"gff file written to {gff}",
    ]

    [log.record(i) for i in info]

    return gff
