#!/usr/bin/env python3

"""
Constants class for postoga. 

postoga is a tool that automates the post-processing of TOGA results.
At its core, this tool takes a TOGA results directory and produces
a series of steps to reduce the amount of manual work required to pre-process
TOGA results for downstream analysis.
"""


__author__ = "Alejandro Gonzales-Irribarren"
__version__ = "0.1.0"
__email__ = "jose.gonzalesdezavala1@unmsm.edu.pe"
__github__ = "https://github.com/alejandrogzi"
__credits__ = ["Bogdan Kirilenko"]


import os


class Constants:
    DESCRIPTION = "postoga is a tool that automates the post-processing of TOGA results."

    TEMP = "temp"

    class ToolNames:
        BED2GTF = "bed2gtf"
        BED2GFF = "bed2gff"

    class FileNames():
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
        SCORES = os.path.join("temp","orthology_scores.tsv")
