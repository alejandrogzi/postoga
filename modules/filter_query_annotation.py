#!/usr/bin/env python3


"""A module to filter the original .bed file based on the query table."""

import os
from typing import Dict, List, Optional, Tuple, Union

import pandas as pd
import polars as pl

from constants import Constants
from logger import Log
from modules.utils import bed_reader

__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "jose.gonzalesdezavala1@unmsm.edu.pe"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.9.3-devel"


def filter_bed(
    togadir: Union[str, os.PathLike],
    outdir: Union[str, os.PathLike],
    table: Union[pd.DataFrame, pl.DataFrame],
    by_class: Optional[Union[str, os.PathLike]],
    by_rel: Optional[Union[str, os.PathLike]],
    threshold: Optional[float],
    paralog: Optional[float],
    bed_path: Union[str, os.PathLike],
    engine: Union[str, os.PathLike] = "pandas",
) -> Tuple:
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

    log = Log.connect(outdir, Constants.FileNames.LOG)
    filtered_bed_path = os.path.join(outdir, Constants.FileNames.FILTERED_BED)
    initial = len(table)

    if threshold:
        if engine != "polars":
            table = table[table["orthology_probability"] >= float(threshold)]
        else:
            table = table.filter(pl.col("orthology_probability") >= threshold)
        log.record(
            f"discarded {initial - len(table)} projections with orthology orthology_probabilitys <{threshold}"
        )

    if by_class:
        edge = len(table)
        if engine != "polars":
            table = table[table["status"].isin(by_class.split(","))]
        else:
            table = table.filter(pl.col("status").is_in(by_class.split(",")))
        log.record(
            f"discarded {edge - len(table)} projections with classes other than {by_class}"
        )

    if by_rel:
        edge = len(table)
        if engine != "polars":
            table = table[table["orthology_relation"].isin(by_rel.split(","))]
        else:
            table = table.filter(pl.col("orthology_relation").is_in(by_rel.split(",")))
        log.record(
            f"discarded {edge - len(table)} projections with orthology_relationships other than {by_rel}"
        )
    if paralog:
        edge = len(table)
        if engine != "polars":
            table = table.groupby("reference_transcript").filter(
                lambda x: (x["orthology_probability"] > float(paralog)).sum() <= 1
            )
        else:
            table = table.filter(
                pl.col("reference_transcript").is_in(
                    table.group_by("reference_transcript")
                    .agg(
                        [
                            (pl.col("orthology_probability") > paralog)
                            .sum()
                            .alias("count_above_paralog")
                        ]
                    )
                    .filter(pl.col("count_above_paralog") <= 1)["helper"]
                )
            )

        log.record(
            f"discarded {edge - len(table)} transcripts with more than 1 chain with orthology probs >{paralog}"
        )

    # read the original .bed file and filter it based on the transcripts table
    if engine != "polars":
        bed = pd.read_csv(bed_path, sep="\t", header=None)
        bed = bed[bed[3].isin(table["projection"])]
        custom_table = table[table["projection"].isin(bed[3])]

        bed.to_csv(
            filtered_bed_path,
            sep="\t",
            header=None,
            index=False,
        )

        stats = [
            custom_table["status"].value_counts().to_dict(),
            custom_table["orthology_relation"].value_counts().to_dict(),
        ]
    else:
        bed = pl.read_csv(
            bed_path,
            separator="\t",
            has_header=False,
        )
        bed = bed.filter(pl.col("column_4").is_in(table["projection"]))
        custom_table = table.filter(pl.col("projection").is_in(bed["projection"]))

        bed.write_csv(filtered_bed_path, include_header=False, separator="\t")

        stats = [
            dict(
                zip(
                    *table.group_by("status")
                    .agg(pl.len())
                    .to_dict(as_series=False)
                    .values()
                )
            ),
            dict(
                zip(
                    *table.group_by("orthology_relation")
                    .agg(pl.len())
                    .to_dict(as_series=False)
                    .values()
                )
            ),
        ]

    info = [
        f"kept {len(bed)} projections after filters, discarded {initial - len(bed)}.",
        f"{len(bed)} projections are coming from {len(custom_table['reference_transcript'].unique())} unique transcripts and {len(custom_table['reference_gene'].unique())} genes",
        f"class stats of new bed: {custom_table['status'].value_counts().to_dict()}",
        f"orthology_relation stats of new bed: {custom_table['orthology_relation'].value_counts().to_dict()}",
        f"filtered bed file written to {f}",
    ]

    [log.record(i) for i in info]

    return (
        filtered_bed_path,
        stats,
        len(custom_table["reference_gene"].unique()),
        custom_table,
    )


def get_stats_from_bed(
    bed_path: str,
    table: Union[pd.DataFrame, pl.DataFrame],
    engine: Union[str, os.PathLike] = "pandas",
) -> Tuple[List[Dict], int]:
    """
    Get the stats of a given bed file

    Parameters
    ----------
        bed_path : str
            The path to the bed file.
        table : pd.DataFrame | pl.DataFrame
            The table to get the stats from.
        engine : str
            The engine to use. Default is "pandas".

    Returns
    -------
        tuple
            The stats and the number of transcripts in the bed file.

    Example
    -------
        >>> from modules.filter_query_annotation import get_stats_from_bed
        >>> get_stats_from_bed("path/to/bed", table)
    """
    bed = bed_reader(bed_path)

    if engine != "polars":
        bed_table = table[table["projection"].isin(bed[3])]
        stats = [
            bed_table["status"].value_counts().to_dict(),
            bed_table["orthology_relation"].value_counts().to_dict(),
        ]
    else:
        bed_table = table.filter(pl.col("projection").is_in(bed["column_4"]))
        stats = [
            dict(
                zip(
                    *table.group_by("status")
                    .agg(pl.len())
                    .to_dict(as_series=False)
                    .values()
                )
            ),
            dict(
                zip(
                    *table.group_by("orthology_relation")
                    .agg(pl.len())
                    .to_dict(as_series=False)
                    .values()
                )
            ),
        ]

    return stats, len(bed_table["reference_gene"].unique())


def unfragment_projections(
    table: Union[pd.DataFrame, pl.DataFrame], togadir: Union[os.PathLike, str], outdir: Union[os.PathLike, str], bed: Union[str, os.PathLike]
) -> Tuple[pd.DataFrame, str]:
    """
    Handles fragmented projections within a BED file by appending unique suffixes
    to duplicate projection IDs and updating a given DataFrame with fragment counts.

    If a projection ID appears multiple times in the BED file (indicating fragments),
    this function modifies the BED file by appending '#FGX' (where X is the fragment
    number) to each instance of the fragmented projection. It also adds a 'fragments'
    column to the input DataFrame, indicating how many times each projection was
    fragmented.

    Parameters
    ----------
    table : pd.DataFrame
        The input DataFrame containing a 'projection' column, which will be
        updated with fragment counts.
    togadir : Union[os.PathLike, str]
        The directory where the original BED file is located and where the
        (potentially) fragmented BED file will be saved.

    Returns
    -------
    tuple
        A tuple containing:
        - pd.DataFrame: The updated DataFrame with the 'fragments' column.
        - str: The path to the (potentially) new fragmented BED file.

    Example
    -------
    >>> import pandas as pd
    >>> import os
    >>> # Assume 'Constants' and 'append_fragment_suffix' are defined
    >>> # Create a dummy BED file
    >>> dummy_bed_content = "chr1\\t10\\t20\\tprojA\\t0\\t+\\n" \\
    ...                     "chr1\\t25\\t35\\tprojB\\t0\\t+\\n" \\
    ...                     "chr1\\t40\\t50\\tprojA\\t0\\t+\\n" \\
    ...                     "chr1\\t55\\t65\\tprojC\\t0\\t+\\n"
    >>> dummy_togadir = "temp_data"
    >>> os.makedirs(dummy_togadir, exist_ok=True)
    >>> with open(os.path.join(dummy_togadir, Constants.FileNames.BED), "w") as f:
    ...     f.write(dummy_bed_content)
    >>>
    >>> # Create a dummy table
    >>> dummy_table = pd.DataFrame({'projection': ['projA', 'projB', 'projC', 'projD'],
    ...                             'value': [10, 20, 30, 40]})
    >>>
    >>> updated_table, bed_path = unfragment_projections(dummy_table, dummy_togadir)
    >>> print(updated_table)
    # Expected output (might vary slightly based on full context of Constants):
    #   projection  value  fragments
    # 0      projA     10          2
    # 1      projB     20          0
    # 2      projC     30          0
    # 3      projD     40          0
    >>> print(bed_path)
    # Expected output: temp_data/fragmented_bed_file.bed
    >>> # Clean up dummy files
    >>> os.remove(os.path.join(dummy_togadir, Constants.FileNames.BED))
    >>> os.remove(os.path.join(dummy_togadir, Constants.FileNames.FRAGMENTED_BED))
    >>> os.rmdir(dummy_togadir)
    """
    bed_path = os.path.join(togadir, bed)
    bed = pd.read_csv(bed_path, sep="\t", header=None)

    projection_count = bed.iloc[:, 3].value_counts()
    fragments = projection_count[projection_count > 1].index

    # INFO: if there are no fragments, does not make sense to
    # update table nor bed
    if len(fragments) > 0:
        fragment_count = {fg: 0 for fg in fragments}
        bed.iloc[:, 3] = bed.apply(
            lambda row: append_fragment_suffix(row, fragment_count), axis=1
        )

        # INFO: add fragment_count as another column in table called "fragments"
        table["fragments"] = [
            fragment_count[projection] if projection in fragment_count.keys() else 0
            for projection in table["projection"]
        ]

        bed.to_csv(
            os.path.join(togadir, Constants.FileNames.FRAGMENTED_BED),
            sep="\t",
            header=None,
            index=False,
        )

        bed_path = os.path.join(togadir, Constants.FileNames.FRAGMENTED_BED)

    table.to_csv(
        os.path.join(outdir, Constants.FileNames.TOGA_TABLE), sep="\t", index=False),
        sep="\t",
        index=False,
    )

    return table, bed_path


def append_fragment_suffix(row, counts):
    """
    Appends a fragment suffix to a projection ID if it is a known fragment.

    This helper function is used within `unfragment_projections` to modify
    projection IDs in the BED file. If a projection ID (from `row[3]`) is
    found in the `counts` dictionary (meaning it's a fragmented projection),
    a suffix like '#FG1', '#FG2', etc., is appended to make it unique.
    The `counts` dictionary is updated to track the next fragment number.

    Parameters
    ----------
    row : pandas.Series
        A row from a pandas DataFrame (typically a BED file row), where
        `row[3]` is expected to be the projection ID.
    counts : dict
        A dictionary where keys are fragmented projection IDs and values are
        the current count of how many times that fragment has been encountered.
        This dictionary is modified in place.

    Returns
    -------
    str
        The modified projection ID with the fragment suffix, or the original
        projection ID if it's not a fragmented one.

    Example
    -------
    >>> # This function is typically used within a pandas .apply() method
    >>> # Example of how it might be used internally:
    >>> import pandas as pd
    >>> counts_dict = {'projA': 0, 'projB': 0}
    >>>
    >>> # First occurrence of projA
    >>> row1 = pd.Series([None, None, None, 'projA', None, None])
    >>> new_id1 = append_fragment_suffix(row1, counts_dict)
    >>> print(f"New ID 1: {new_id1}, Counts: {counts_dict}")
    # Expected: New ID 1: projA#FG1, Counts: {'projA': 1, 'projB': 0}
    >>>
    >>> # Second occurrence of projA
    >>> row2 = pd.Series([None, None, None, 'projA', None, None])
    >>> new_id2 = append_fragment_suffix(row2, counts_dict)
    >>> print(f"New ID 2: {new_id2}, Counts: {counts_dict}")
    # Expected: New ID 2: projA#FG2, Counts: {'projA': 2, 'projB': 0}
    >>>
    >>> # A non-fragmented ID
    >>> row3 = pd.Series([None, None, None, 'projX', None, None])
    >>> new_id3 = append_fragment_suffix(row3, counts_dict)
    >>> print(f"New ID 3: {new_id3}, Counts: {counts_dict}")
    # Expected: New ID 3: projX, Counts: {'projA': 2, 'projB': 0}
    """
    val = row[3]
    if val in counts:
        count = counts[val]
        idx = counts[val] = count + 1
        return f"{val}#FG{idx}"
    else:
        return val
