"""
This module processes orthology classification, loss summary, and scoring data
to create a unified table of gene orthology relationships.
"""

import logging
import os
from typing import Dict, List, Optional, Tuple, Union

import pandas as pd

__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "alejandrxgzi@gmail.com"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.12"

ORTHOLOGY_CLASSIFICATION_COLS = [
    "reference_gene",
    "reference_transcript",
    "query_gene",
    "query_transcript",
    "orthology_class",
]

LOSS_SUMMARY_COLS = [
    "level",
    "query_transcript",
    "loss_status",
]

ORTHOLOGY_SCORE_COLS = [
    "transcript",
    "chain",
    "orthology_score",
]

BED12_COLS = [
    "chr",
    "start",
    "end",
    "id",
    "score",
    "strand",
    "cds_start",
    "cds_end",
    "rgb",
    "exon_count",
    "exon_sizes",
    "exon_starts",
]

OUTPUT_TABLE_COLS = [
    "reference_gene",
    "reference_transcript",
    "query_gene",
    "query_transcript",
    "orthology_class",
    "loss_status",
    "orthology_score",
]

PathLike = Union[str, os.PathLike]
LOGGER = logging.getLogger("postoga")


def load_tsv_with_columns(
    filepath: PathLike, column_names: List[str], skip_header: bool = True
) -> pd.DataFrame:
    """
    Load a TSV file with specified column names.

    Args:
        filepath: Path to the TSV file
        column_names: List of column names to assign
        skip_header: Whether to skip the first row

    Returns:
        DataFrame with the loaded data

    Example:
        >>> df = load_tsv_with_columns("data.tsv", ["col1", "col2"])
    """
    return pd.read_csv(
        filepath, sep="\t", names=column_names, skiprows=1 if skip_header else 0
    )


def load_orthology_classification(filepath: PathLike) -> pd.DataFrame:
    """
    Load orthology classification data.

    Args:
        filepath: Path to orthology_classification.tsv

    Returns:
        DataFrame with orthology classifications

    Example:
        >>> df = load_orthology_classification("orthology_classification.tsv")
    """
    return load_tsv_with_columns(filepath, ORTHOLOGY_CLASSIFICATION_COLS)


def load_and_process_scores(filepath: PathLike) -> pd.DataFrame:
    """
    Load and process orthology scores, creating query_transcript identifiers.

    Args:
        filepath: Path to orthology_scores.tsv

    Returns:
        DataFrame with query_transcript, orthology_score, and transcript columns

    Example:
        >>> df = load_and_process_scores("orthology_scores.tsv")
    """
    scores = load_tsv_with_columns(filepath, ORTHOLOGY_SCORE_COLS)
    scores["query_transcript"] = (
        scores["transcript"] + "#" + scores["chain"].astype(str)
    )
    return scores[["query_transcript", "orthology_score", "transcript"]]


def load_and_filter_loss_summary(
    filepath: PathLike, level_filter: str = "PROJECTION"
) -> pd.DataFrame:
    """
    Load loss summary data and filter by level.

    Args:
        filepath: Path to loss_summary.tsv
        level_filter: Level to filter by (default: "PROJECTION")

    Returns:
        DataFrame with filtered loss summary and r_transcript column

    Example:
        >>> df = load_and_filter_loss_summary("loss_summary.tsv")
    """
    loss = load_tsv_with_columns(filepath, LOSS_SUMMARY_COLS)
    loss = loss.query(f"level == '{level_filter}'").copy()
    loss["r_transcript"] = loss["query_transcript"].str.split("#").str[:2].str.join("#")
    return loss


def load_query_genes_mapping(filepath: PathLike) -> Dict[str, str]:
    """
    Load query genes mapping from projection to query_gene.

    Args:
        filepath: Path to query_genes.tsv

    Returns:
        Dictionary mapping projection IDs to query gene names

    Example:
        >>> mapping = load_query_genes_mapping("query_genes.tsv")
    """
    return (
        pd.read_csv(filepath, sep="\t").set_index("projection")["query_gene"].to_dict()
    )


def load_query_annotation(filepath: PathLike) -> pd.DataFrame:
    """
    Load query annotation BED file and extract helper transcript IDs.

    Args:
        filepath: Path to query_annotation.bed

    Returns:
        DataFrame with BED12 data and helper column

    Example:
        >>> df = load_query_annotation("query_annotation.bed")
    """
    annotation = pd.read_csv(filepath, sep="\t", header=None, names=BED12_COLS)
    annotation["helper"] = annotation["id"].str.split("$").str[0]
    return annotation


def create_gene_mapping(
    df: pd.DataFrame, key_col: str, value_col: str
) -> Dict[str, str]:
    """
    Create a mapping dictionary from two columns, dropping NaN values.

    Args:
        df: Source DataFrame
        key_col: Column name for dictionary keys
        value_col: Column name for dictionary values

    Returns:
        Dictionary mapping keys to values

    Example:
        >>> mapping = create_gene_mapping(df, 'transcript', 'gene')
    """
    return df[[key_col, value_col]].dropna().set_index(key_col)[value_col].to_dict()


def fill_gene_columns_with_mapping(
    df: pd.DataFrame,
    mapping: Dict[str, str],
    transcript_col: str = "reference_transcript",
) -> pd.DataFrame:
    """
    Fill NaN values in reference_gene and query_gene using a transcript mapping.

    Args:
        df: DataFrame to update
        mapping: Dictionary mapping transcripts to genes
        transcript_col: Column to use for mapping lookup

    Returns:
        DataFrame with filled gene columns

    Example:
        >>> df = fill_gene_columns_with_mapping(df, gene_map)
    """
    df = df.copy()
    df["reference_gene"] = df["reference_gene"].fillna(df[transcript_col].map(mapping))
    df["query_gene"] = df["query_gene"].fillna(df[transcript_col].map(mapping))
    return df


def merge_and_fill_orthology_data(
    orthology: pd.DataFrame, loss: pd.DataFrame, scores: pd.DataFrame
) -> pd.DataFrame:
    """
    Merge orthology, loss, and score data with intelligent filling of missing values.

    Args:
        orthology: Orthology classification DataFrame
        loss: Loss summary DataFrame
        scores: Orthology scores DataFrame

    Returns:
        Merged DataFrame with filled reference and query information

    Example:
        >>> merged = merge_and_fill_orthology_data(orth_df, loss_df, score_df)
    """
    # Initial merge with loss data
    table = pd.merge(orthology, loss, on="query_transcript", how="outer").fillna(
        {"reference_transcript": loss["r_transcript"]}
    )

    # Fill genes using reference_transcript mapping
    gene_mapping = create_gene_mapping(table, "reference_transcript", "reference_gene")
    table = fill_gene_columns_with_mapping(table, gene_mapping)

    # Merge with scores
    table = table.merge(scores, on="query_transcript", how="outer")
    table = table.fillna({"reference_transcript": table["transcript"]})

    # Fill genes again with updated mapping
    gene_mapping = create_gene_mapping(table, "reference_transcript", "reference_gene")
    table = fill_gene_columns_with_mapping(table, gene_mapping)

    return table


def apply_query_gene_overrides(
    df: pd.DataFrame, query_genes_mapping: Dict[str, str]
) -> pd.DataFrame:
    """
    Override query_gene values where query_transcript exists in mapping.

    Args:
        df: DataFrame with query_transcript and query_gene columns
        query_genes_mapping: Dictionary mapping transcripts to gene names

    Returns:
        DataFrame with overridden query_gene values

    Example:
        >>> df = apply_query_gene_overrides(df, gene_mapping)
    """
    df = df.copy()
    mask = df["query_transcript"].isin(query_genes_mapping.keys())
    df.loc[mask, "query_gene"] = df.loc[mask, "query_transcript"].map(
        query_genes_mapping
    )
    return df


def extract_isoform_mappings(
    annotation: pd.DataFrame, orthology_table: pd.DataFrame
) -> pd.DataFrame:
    """
    Extract isoform to query_gene mappings from annotation and orthology data.

    Args:
        annotation: Query annotation DataFrame with 'id' and 'helper' columns
        orthology_table: Orthology table with query_gene and query_transcript

    Returns:
        DataFrame mapping isoform IDs to query genes

    Example:
        >>> isoforms = extract_isoform_mappings(annot_df, orth_df)
    """
    return annotation[["id", "helper"]].merge(
        orthology_table[["query_gene", "query_transcript"]],
        left_on="helper",
        right_on="query_transcript",
    )[["id", "query_gene"]]


def table_builder(
    query_annotation: Union[str, PathLike],
    loss_summary: Union[str, PathLike],
    orthology_classification: Union[str, PathLike],
    orthology_scores: Union[str, PathLike],
    query_genes: Union[str, PathLike],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Execute the complete orthology data processing pipeline.

    Args:
        toga_dir: Base directory containing all input files

    Returns:
        Tuple of (orthology_table, isoform_mappings)
        - orthology_table: Complete merged orthology data
        - isoform_mappings: Isoform ID to query gene mappings

    Example:
        >>> table, isoforms = process_orthology_pipeline("./data")
    """
    # Load all data sources
    orthology = load_orthology_classification(orthology_classification)
    scores = load_and_process_scores(orthology_scores)
    loss = load_and_filter_loss_summary(loss_summary)
    query_genes_mapping = load_query_genes_mapping(query_genes)
    query_annotation = load_query_annotation(query_annotation)

    # Process and merge orthology data
    orthology_table = merge_and_fill_orthology_data(orthology, loss, scores)
    orthology_table = apply_query_gene_overrides(orthology_table, query_genes_mapping)
    orthology_table = orthology_table[OUTPUT_TABLE_COLS]

    # Fill and change query_genes only for fragmented projections
    orthology_table = fill_query_genes_for_fragmented_projections(orthology_table)

    # Extract isoform mappings
    isoform_mappings = extract_isoform_mappings(query_annotation, orthology_table)

    return orthology_table, isoform_mappings


def fill_query_genes_for_fragmented_projections(table: pd.DataFrame) -> pd.DataFrame:
    """
    Fill query_genes for fragmented projections.

    Args:
        table: Orthology table with query_gene and query_transcript columns

    Returns:
        DataFrame with filled query_genes

    Example:
        >>> table = fill_query_genes_for_fragmented_projections(table)
    """
    table = table.copy()

    # Create mask for rows with "$"
    mask = table["query_transcript"].str.contains("$", regex=False)

    # Process only rows that need updating
    if mask.any():
        fragments = table.loc[mask, "query_transcript"].str.split("$").str[-1]
        table.loc[mask, "query_gene"] = table.loc[mask, "query_gene"] + "_" + table.loc[mask, "query_transcript"]

    return table

def filter_query_annotation(
    table: pd.DataFrame,
    by_orthology_class: Optional[Union[str, os.PathLike]],
    by_loss_status: Optional[Union[str, os.PathLike]],
    by_orthology_score: Optional[float],
    by_paralog_score: Optional[float],
    query_annotation: Union[str, os.PathLike],
    outdir: Union[os.PathLike[str], str],
) -> Tuple[Union[str, os.PathLike], pd.DataFrame]:
    """
    Filters the original .bed file to produce a custom filtered file

    Parameters
    ----------
        togadir : Union[str, os.PathLike]|Union[str, os.PathLike]
            The path to the TOGA results directory.
        outdir : Union[str, os.PathLike]|Union[str, os.PathLike]
            The path to the output directory.
        table : pd.DataFrame
            The query table.
        by_class : list
            The classes to filter the table by.
        by_rel : list
            The orthology_relationships to filter the table by.
        threshold : float
            The orthology orthology_probability threshold.
        paralog : float
            The paralogy orthology_probability threshold.

    Returns
    -------
        tuple
            The filtered table and the custom table.

    Example
    -------
        >>> from modules.filter_query_annotation import filter_bed
    """

    tmp = query_annotation.replace(".bed", ".filtered.bed")
    initial_table_size = len(table)
    LOGGER.debug("Starting BED filtering with %d projections", initial_table_size)

    if by_orthology_score is not None:
        edge = len(table)
        table = table[table["orthology_score"] >= float(by_orthology_score)]
        discarded = edge - len(table)
        if discarded:
            LOGGER.debug(
                "Discarded %d projections with orthology_score <%s",
                discarded,
                by_orthology_score,
            )

    if by_orthology_class:
        edge = len(table)
        allowed_classes = [
            klass.strip() for klass in by_orthology_class.split(",") if klass.strip()
        ]
        table = table[table["orthology_class"].isin(allowed_classes)]
        discarded = edge - len(table)
        if discarded:
            LOGGER.debug(
                "Discarded %d projections outside classes %s",
                discarded,
                ",".join(allowed_classes),
            )

    if by_loss_status:
        edge = len(table)
        allowed_relationships = [
            relation.strip()
            for relation in by_loss_status.split(",")
            if relation.strip()
        ]
        table = table[table["loss_status"].isin(allowed_relationships)]
        discarded = edge - len(table)
        if discarded:
            LOGGER.debug(
                "Discarded %d projections outside relationships %s",
                discarded,
                ",".join(allowed_relationships),
            )

    if by_paralog_score is not None:
        edge = len(table)
        if "orthology_score" not in table.columns:
            LOGGER.warning(
                "Paralog filtering skipped: 'orthology_score' column not found."
            )
        else:
            table = table.groupby("reference_transcript").filter(
                lambda x: (x["orthology_score"] > float(by_paralog_score)).sum() <= 1
            )
            discarded = edge - len(table)
            if discarded:
                LOGGER.debug(
                    "Discarded %d transcripts with paralog score >%s",
                    discarded,
                    by_paralog_score,
                )

    # read the original .bed file and filter it based on the transcripts table
    query_annotation_df = load_query_annotation(query_annotation)
    filtered_query_annotation = query_annotation_df[
        query_annotation_df["helper"].isin(table["query_transcript"])
    ]
    filtered_table = table[
        table["query_transcript"].isin(filtered_query_annotation["helper"])
    ]

    kept = len(filtered_query_annotation)
    discarded = initial_table_size - kept
    unique_transcripts = (
        filtered_table["reference_transcript"].nunique()
        if "reference_transcript" in filtered_table
        else 0
    )
    unique_genes = (
        filtered_table["reference_gene"].nunique()
        if "reference_gene" in filtered_table
        else 0
    )
    class_stats = (
        filtered_table["orthology_class"].value_counts().to_dict()
        if "orthology_class" in filtered_table
        else {}
    )
    status_stats = (
        filtered_table["loss_status"].value_counts().to_dict()
        if "loss_status" in filtered_table
        else {}
    )

    filtered_query_annotation[BED12_COLS].to_csv(
        tmp, sep="\t", index=False, header=False
    )

    info_messages = [
        f"kept {kept} projections after filters, discarded {discarded}.",
        f"""{kept} projections from {unique_transcripts} unique tx and {
            unique_genes
        } genes""",
        f"class stats of new bed: {class_stats}",
        f"orthology_relation stats of new bed: {status_stats}",
        f"filtered bed file written to {tmp}",
    ]
    for message in info_messages:
        LOGGER.info(message)

    return (tmp, filtered_table)
