#!/usr/bin/env python3


"""
Constants class for postoga. 

postoga is a tool that automates the post-processing of TOGA results.
At its core, this tool takes a TOGA results directory and produces
a series of steps to reduce the amount of manual work required to pre-process
TOGA results for downstream analysis.
"""


import os


__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "jose.gonzalesdezavala1@unmsm.edu.pe"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.4.0-devel"
__credits__ = ["Bogdan Kirilenko"]


class Constants:
    DESCRIPTION = (
        "postoga is a tool that automates the post-processing of TOGA results."
    )

    TEMP = "temp"
    ORTHOLOGY_TYPE = {
        "one2one": "o2o",
        "one2many": "o2m",
        "many2one": "m2o",
        "many2many": "m2m",
        "one2zero": "o2z",
    }
    ORDER = {"I": 1, "PI": 2, "UL": 3, "L": 4, "M": 5, "PM": 6, "PG": 7, "abs": 8}

    class ToolNames:
        BED2GTF = "bed2gtf"
        BED2GFF = "bed2gff"
        COMPLEASM = "compleasm"

    class FileNames:
        ORTHOLOGY = "orthology_classification.tsv"
        BED = "query_annotation.bed"
        CODON = "codon.fasta"
        PROTEIN = "prot.fasta"
        ISOFORMS = os.path.join("temp", "isoforms.tsv")
        OWNED_ISOFORMS = "isoforms.txt"
        FILTERED_BED = "filtered.bed"
        GTF = f"{BED.split('.')[0]}.gtf"
        GFF = f"{BED.split('.')[0]}.gff"
        FILTERED_GTF = f"{FILTERED_BED.split('.')[0]}.gtf"
        FILTERED_GFF = f"{FILTERED_BED.split('.')[0]}.gff"
        NUCLEOTIDE = "nucleotide.fasta"
        CLASS = "loss_summ_data.tsv"
        SCORES = os.path.join("temp", "orthology_scores.tsv")
        LOG = "postoga.log"
        QUALITY = os.path.join("temp", "transcript_quality.tsv")
        ANCESTRAL = os.path.join("./supply", "Ancestral_placental.txt")

    class Commands:
        COMMIT = "git rev-parse --short HEAD"
        BRANCH = "git rev-parse --abbrev-ref HEAD"

    class Metadata:
        BED2GTF_METADATA = "https://github.com/alejandrogzi/bed2gtf"
        BED2GFF_METADATA = "https://github.com/alejandrogzi/bed2gff"
