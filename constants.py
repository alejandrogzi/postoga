#!/usr/bin/env python3

__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "jose.gonzalesdezavala1@unmsm.edu.pe"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.10.1"
__credits__ = ["Bogdan Kirilenko", "Yury Malovichko"]


class Constants:
    DESCRIPTION = (
        "postoga is a tool that automates the post-processing of TOGA results."
    )

    TEMP = "temp"

    class ToolNames:
        BED2GTF = "bed2gtf"
        BED2GFF = "bed2gff"

    class FileNames:
        ORTHOLOGY = "orthology_classification.tsv"
        LOSS = "loss_summary.tsv"
        SCORES = "orthology_scores.tsv"
        INACTIVATING_MUTATIONS = "inactivating_mutations.tsv"
        PARALOGS = "meta/paralogous_projections.tsv"
        CODON = "codon_aln.fa"
        FILTERED_CODON = "codon_aln.filtered.fa"
        PROTEIN = "protein_aln.fa"
        QUERY_GENES = "query_genes.tsv"
        FILTERED_PROTEIN = "protein_aln.filtered.fa"
        EXONS = "exon_aln.fa"
        OWNED_ISOFORMS = "postoga_isoforms.txt"
        NUCLEOTIDE = "nucleotide.fasta"  # [DEPRECATED IN TOGA2]
        LOG = "postoga.log"
        TOGA_TABLE = "toga.table.gz"
        BED = "query_annotation.bed"
        BED_UTR = "query_annotation.with_utrs.bed"
        FRAGMENTED_BED = f"{BED.rsplit('.', 1)[0]}.fragmented.bed"
        FILTERED_BED = f"{BED.rsplit('.', 1)[0]}.filtered.bed"
        GTF = f"{BED.rsplit('.', 1)[0]}.gtf"
        GFF = f"{BED.rsplit('.', 1)[0]}.gff"
        FILTERED_GTF = f"{FILTERED_BED.split('.')[0]}.gtf"
        FILTERED_GFF = f"{FILTERED_BED.split('.')[0]}.gff"
        POSTOGA_OUTPUT_DIRECTORY = "POSTOGA"

    class Commands:
        COMMIT = "git rev-parse --short HEAD"
        BRANCH = "git rev-parse --abbrev-ref HEAD"

    class Metadata:
        BED2GTF_METADATA = "https://github.com/alejandrogzi/bed2gtf"
        BED2GFF_METADATA = "https://github.com/alejandrogzi/bed2gff"
