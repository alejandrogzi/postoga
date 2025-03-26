#!/usr/bin/env python3

import os

import supply

try:
    import importlib.resources as resources
except ImportError:
    import importlib_resources as resources

__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "jose.gonzalesdezavala1@unmsm.edu.pe"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.9.3-devel"
__credits__ = ["Bogdan Kirilenko", "Yury Malovichko"]


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
        # "high_confidence": "#0073ae",
        # "partial": "#6ab1e0",
        # "average_confidence": "#009bd9",
        # "low_confidence": "#cf5058",
    }
    CATEGORY_ORDER = ["I", "PI", "UL", "L", "M", "PM", "o2o", "m2m", "m2o", "o2m"]
    STACKED_COLUMN_NAMES = {
        "Count_0": "Orthology\nclass",
        "Count_1": "Orthology\nrelationship",
        "Count_2": "Orthology\nconfidence",
    }
    ANCESTRAL_CATEGORY = {"mut": ["UL", "L"], "missing": ["PI", "M", "PM", "PG", "NF"]}
    ANCESTRAL_NGENES_ENSEMBL = 18430
    ANCESTRAL_NGENES_ENTREZ = 18268
    ANCESTRAL_NGENES_NAME = 18359
    ANCESTRAL_NGENES = {
        "ensembl": ANCESTRAL_NGENES_ENSEMBL,
        "entrez": ANCESTRAL_NGENES_ENTREZ,
        "gene_name": ANCESTRAL_NGENES_NAME,
    }
    SUPERORDER = ["Laurasiatheria", "Euarchontoglires", "User"]
    SUPERORDER_COLORS = {
        "Laurasiatheria": "black",
        "Euarchontoglires": "blue",
        "Other": "grey",
        "User": "red",
    }
    SPECIES_DEFAULT = "human"
    SRC_DEFAULT = "ensembl"
    PLOTSTAMP = """Generated on {} by postoga \nversion: {}, branch: {}, commit: {}.\n
    This report provides a basic analysis of the data and results
    obtained by TOGA and is intended to be used as a preliminary
    analisys step. For more information, updates, bugs, suggestions
    or any other inquire, please visit our GitHub repository at
    github.com/alejandrogzi/postoga."""
    HIST_NBINS = 50
    PHYLO_DEFAULT = "mammals"
    BUSCO_DBS_BASE = {"eukaryota": "eukaryota_odb10", "vertebrata": "vertebrata_odb10"}
    BUSCO_DBS_MAMMALS = {
        "mammalia": "mammalia_odb10",
        "eutheria": "eutheria_odb10",
    }  # Carnivora
    BUSCO_DBS_BIRDS = {"aves": "aves_odb10"}  # Passeriformes

    class ToolNames:
        BED2GTF = "bed2gtf"
        BED2GFF = "bed2gff"
        NOEL = "noel"

    class DirNames:
        FIGURES = "POSTOGA_FIGURES"

    class FigNames:
        PLOTTING_FORMAT = "jpeg"  # choices = ["pdf", "png", "svg", "jpg", "jpeg", "tif", "tiff", "eps", "ps"]
        LENGTHS_PLOT = f"ortholog_lengths.{PLOTTING_FORMAT}"
        BARPLOT = f"qual_barplot.{PLOTTING_FORMAT}"
        SCORE_PLOT = f"orthology_scores.{PLOTTING_FORMAT}"
        ANNOTATION_BOXPLOT = f"annotation_boxplot.{PLOTTING_FORMAT}"
        ANCESTRAL_BARPLOT = f"ancestral_barplot.{PLOTTING_FORMAT}"
        ANCESTRAL_SCATTER = f"ancestral_scatterplot.{PLOTTING_FORMAT}"
        BUSCO_BARPLOT = f"busco_completeness.{PLOTTING_FORMAT}"

    class FileNames:
        SUPPLY_FOLDER = resources.files(supply)
        ORTHOLOGY = "orthology_classification.tsv"
        LOSS = "loss_summary.tsv"
        SCORES = "orthology_scores.tsv"
        INACTIVATING_MUTATIONS = "inactivating_mutations.tsv"
        PARALOGS = "meta/paralogous_projections.tsv"
        BED = "query_annotation.bed"
        CODON = "codon_aln.fa"
        FILTERED_CODON = "codon_aln.filtered.fa"
        PROTEIN = "protein_aln.fa"
        QUERY_GENES = "query_genes.tsv"
        FILTERED_PROTEIN = "protein_aln.filtered.fa"
        EXONS = "exon_aln.fa"
        OWNED_ISOFORMS = "postoga_isoforms.txt"
        FILTERED_BED = "filtered.bed"
        GTF = f"{BED.split('.')[0]}.gtf"
        GFF = f"{BED.split('.')[0]}.gff"
        FILTERED_GTF = f"{FILTERED_BED.split('.')[0]}.gtf"
        FILTERED_GFF = f"{FILTERED_BED.split('.')[0]}.gff"
        NUCLEOTIDE = "nucleotide.fasta"  # [DEPRECATED IN TOGA2]
        LOG = "postoga.log"
        QUALITY = os.path.join("temp", "transcript_quality.tsv")
        ANCESTRAL = SUPPLY_FOLDER.joinpath("Ancestral_placental_complete.txt")
        HAPLOTYPE = "merged_assemblies.txt"
        MAMMALS = SUPPLY_FOLDER.joinpath("mammal_genes_template.txt")
        BIRDS = SUPPLY_FOLDER.joinpath("birds_genes_template.txt")
        LOGO_IMG = SUPPLY_FOLDER.joinpath("postoga_logo.png")
        FONT = SUPPLY_FOLDER.joinpath("font/Arial.ttf")
        PDF = "POSTOGA_REPORT.pdf"
        LENGTHS = "ortholog_lengths.txt"
        TOGA_TABLE = "toga.table"

    class Commands:
        COMMIT = "git rev-parse --short HEAD"
        BRANCH = "git rev-parse --abbrev-ref HEAD"

    class Metadata:
        BED2GTF_METADATA = "https://github.com/alejandrogzi/bed2gtf"
        BED2GFF_METADATA = "https://github.com/alejandrogzi/bed2gff"
        NOEL_METADATA = "https://github.com/alejandrogzi/noel"
