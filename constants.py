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
__version__ = "0.6.0-devel"
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
    ORDER = {
        "I": 1,
        "PI": 2,
        "UL": 3,
        "L": 4,
        "M": 5,
        "PM": 6,
        "PG": 7,
        "N": 8,
        "NF": 9,
    }
    TAXA_ORDER = [
        "Hominoidea",
        "Scandentia",
        "Dermoptera",
        "Rodentia",
        "Lagomorpha",
        "Tylopoda",
        "Whippomorpha",
        "Ruminantia",
        "Suina",
        "Perissodactyla",
        "Carnivora",
        "Pholidota",
        "Chiroptera",
        "Eulipotyphla",
        "Afrotheria",
        "Xenarthra",
        "Marsupialia",
        "Monotremata",
        "Galloanserae",
        "Columbiformes",
        "Cuculimorphae",
        "Caprimulgimorphae",
        "Gruiformes",
        "Charadriiformes",
        "Aesquomithes",
        "Strigiformes",
        "Piciformes",
        "Coraciiformes",
        "Accipitriformes",
        "Falconiformes",
        "Psittaciformes",
        "Passeriformes",
        "Palaeognathae",
    ]
    CATEGORY_COLORS = {
        "I": "#0073ae",
        "UL": "#6ab1e0",
        "L": "#89171a",
        "PM": "#cf5058",
        "M": "#ec212a",
        "PI": "#009bd9",
        "PG": "black",
        "o2o": "#0073ae",
        "m2m": "#009bd9",
        "m2o": "#6ab1e0",
        "o2m": "#cf5058",
        "high_confidence": "#0073ae",
        "partial": "#6ab1e0",
        "average_confidence": "#009bd9",
        "low_confidence": "#cf5058",
    }
    STACKED_COLUMN_NAMES = {
        "Count_0": "Orthology\nclass",
        "Count_1": "Orthology\nrelationship",
        "Count_2": "Orthology\nconfidence",
    }
    ANCESTRAL_CATEGORY = {"mut": ["UL", "L"], "missing": ["PI", "M", "PM", "PG", "NF"]}
    ANCESTRAL_NGENES = 18430
    SUPERORDER = ["Laurasiatheria", "Euarchontoglires", "User"]
    SUPERORDER_COLORS = {
        "Laurasiatheria": "black",
        "Euarchontoglires": "blue",
        "Other": "grey",
        "User": "red",
    }
    SPECIES_DEFAULT = "human"
    PLOTSTAMP = """Generated on {} by postoga \nversion: {}, branch: {}, commit: {}.\n
    This report provides a basic analysis of the data and results 
    obtained by TOGA and is intended to be used as a preliminary 
    analisys step. For more information, updates, bugs, suggestions 
    or any other inquire, please visit our GitHub repository at
    github.com/alejandrogzi/postoga."""

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
        HAPLOTYPE = "merged_assemblies.txt"
        MAMMALS = os.path.join("./supply", "mammal_genes_template.txt")
        BIRDS = os.path.join("./supply", "birds_genes_template.txt")
        PDF = "POSTOGA_REPORT.pdf"

    class Commands:
        COMMIT = "git rev-parse --short HEAD"
        BRANCH = "git rev-parse --abbrev-ref HEAD"

    class Metadata:
        BED2GTF_METADATA = "https://github.com/alejandrogzi/bed2gtf"
        BED2GFF_METADATA = "https://github.com/alejandrogzi/bed2gff"
