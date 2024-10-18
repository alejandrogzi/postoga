#!/usr/bin/env python3

"""Plotter module for postoga."""

import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
from collections import defaultdict
from matplotlib import ticker
import matplotlib.font_manager as font_manager
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
from datetime import datetime
import matplotlib.gridspec as gridspec
from constants import Constants
from logger import Log
from modules.utils import shell


__author__ = "Alejandro Gonzales-Irribarren"
__email__ = "jose.gonzalesdezavala1@unmsm.edu.pe"
__github__ = "https://github.com/alejandrogzi"
__version__ = "0.9.3-devel"


def set_font():
    font_dir = Constants.FileNames.FONT  # "./supply/font/Arial.ttf"
    font_manager.fontManager.addfont(font_dir)
    plt.rcParams["font.family"] = "Arial"


def bar_setup(ax):
    ax.spines["right"].set_color("none")
    ax.spines["top"].set_color("none")
    ax.spines["bottom"].set_color("none")
    ax.spines["left"].set_linewidth(1.5)
    ax.xaxis.set_major_locator(ticker.NullLocator())
    ax.yaxis.set_ticks_position("left")
    ax.tick_params(which="major", width=1.00, direction="in")
    ax.tick_params(which="major", length=7, direction="in")
    ax.tick_params(which="minor", width=0.75, direction="in")
    ax.tick_params(which="minor", length=3.5, direction="in")
    ax.set_ylim(0, 100)
    ax.patch.set_alpha(0.0)
    ax.yaxis.set_major_locator(ticker.MultipleLocator(20))
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(10))

    return ax


def make_boxplot_annotated_genes(
    path: str, reference: str, ngenes: int
):  # reference must be "human", "mouse" or "chicken"; ngenes: number of genes annotated
    set_font()
    fig, ax = plt.subplots(figsize=(9, 4))

    if reference != "chicken":
        df = pd.read_csv(
            Constants.FileNames.MAMMALS, sep="\t", thousands=","
        )  # make it a constant
    else:
        df = pd.read_csv(Constants.FileNames.BIRDS, sep="\t", thousands=",")

    taxonomy = []
    for lineage in df["Taxonomic Lineage"]:
        found_in_orders = False

        lin = [x.strip() for x in lineage.split(";")]

        for order in Constants.TAXA_ORDER:
            if order in lin:
                found_in_orders = (
                    True  # Set the flag to True if any element is found in orders
                )
                taxonomy.append(order)
                break

        if not found_in_orders:
            if "Primates" in lin:
                taxonomy.append("Other Primates")
            elif "Diprotodontia" in lin:
                taxonomy.append("Diprotodontia")
            else:
                taxonomy.append("Other")

    df["taxonomy"] = taxonomy

    if reference != "chicken":
        df.loc[len(df.index), ["Species", "ref_human", "ref_mouse", "taxonomy"]] = (
            "User",
            ngenes,
            ngenes,
            "Your assembly",
        )
        user_place = 21
        if reference == "human":
            ref_genes = 19464
            title_place = 20300
        else:
            ref_genes = 22258
            title_place = 23600
    else:
        df.loc[len(df.index), ["Species", "ref_chicken", "taxonomy"]] = (
            "User",
            ngenes,
            "Your assembly",
        )
        ref_genes = 18039
        user_place = 15
        title_place = 19000

    # Create a dictionary to map taxonomy to ref_human values
    taxonomy_dict = {}
    for taxon, group in df.groupby("taxonomy"):
        taxonomy_dict[taxon] = group[f"ref_{reference}"]

    # Create a list of ref_human values for each taxonomy
    ref_data = [values for values in taxonomy_dict.values()]

    # Create a list of taxonomy labels
    taxonomy_labels = list(taxonomy_dict.keys())

    for i, values in enumerate(ref_data):
        if taxonomy_labels[i] != "Your assembly":
            x = np.random.normal(i + 1, 0.04, size=len(values))
            ax.scatter(
                x, values, label=f"Taxonomy {taxonomy_labels[i]}", s=10, alpha=0.5
            )  # Adjust 's' for point size and 'alpha' for transparency
        else:
            ax.scatter(
                user_place,
                values,
                label=f"Taxonomy {taxonomy_labels[i]}",
                s=13,
                alpha=1,
                c="red",
            )

    ax.boxplot(ref_data, labels=taxonomy_labels, showfliers=False)

    ax.set_xticklabels(taxonomy_labels, rotation=45, ha="right", rotation_mode="anchor")
    ax.legend("", frameon=False)

    for value in [
        8000,
        10000,
        12000,
        14000,
        16000,
        18000,
    ]:  # You can specify the values where you want to add lines
        ax.axhline(
            value,
            xmin=0,
            xmax=len(taxonomy_labels) + 1,
            linewidth=1,
            color="gray",
            linestyle="--",
            dashes=(5, 7),
            alpha=0.3,
        )

    ax.axvline(user_place, color="#fffedd", linewidth=23, alpha=0.5)
    ax.axhline(ref_genes, color="b")
    ax.get_xticklabels()[user_place - 1].set_color("red")
    ax.set_ylabel("No. of annotated\northologous genes")
    ax.text(
        0.5,
        title_place,
        f"{reference} as reference: {ref_genes:,} input genes",
        fontsize=11,
    )

    plt.savefig(
        os.path.join(
            path, Constants.DirNames.FIGURES, Constants.FigNames.ANNOTATION_BOXPLOT
        ),
        format=Constants.FigNames.PLOTTING_FORMAT,
        bbox_inches="tight",
        dpi=300,
    )


def make_query_barplot_stats(path: str, dbs: list):
    fig, (ax1, ax2) = plt.subplots(1, 2)
    axs = [ax1, ax2]

    for idx, db in enumerate(dbs):
        db = [{key: x[key] for key in Constants.CATEGORY_ORDER if key in x} for x in db]
        count = 0
        dfs = []
        for x in db:
            df = pd.DataFrame.from_dict(
                x, orient="index", columns=[f"Count_{count}"]
            ).transpose()
            df = df.div(df.sum(axis=1), axis=0) * 100
            df.fillna(0, inplace=True)
            dfs.append(df)
            count += 1

        combined_df = pd.concat(dfs, axis=1)
        combined_df.fillna(0, inplace=True)
        df = combined_df
        x = df.index
        bottom = None

        for column in df.columns:
            heights = df[column]
            axs[idx].bar(
                x,
                heights,
                bottom=bottom,
                label=column,
                color=Constants.CATEGORY_COLORS[column],
            )

            if bottom is None:
                bottom = heights
            else:
                bottom += heights

        legend_labels = defaultdict(list)

        for category, color in Constants.CATEGORY_COLORS.items():
            if category in df.columns:
                legend_labels[color].append(category)

        custom_legend_labels = [
            " / ".join(categories) for color, categories in legend_labels.items()
        ]

        axs[idx].set_ylabel("Proportion (%)")

        if idx > 0:
            axs[idx].set_title("Filtered annotation", ha="center", fontsize=11)
        else:
            axs[idx].set_title("Base annotation", ha="center", fontsize=11)

            axs[idx].legend(
                custom_legend_labels,
                bbox_to_anchor=(1.85, -0.2),
                frameon=False,
                handlelength=1,
                handleheight=1,
                fontsize="small",
                ncol=4,
            )

        bar_setup(axs[idx])

        for i, col in enumerate(combined_df.transpose().columns):
            axs[idx].text(
                i,
                -10.5,
                f"{Constants.STACKED_COLUMN_NAMES[col]}",
                ha="center",
                va="center",
                fontsize=10,
            )

    plt.savefig(
        os.path.join(path, Constants.DirNames.FIGURES, Constants.FigNames.BARPLOT),
        format=Constants.FigNames.PLOTTING_FORMAT,
        bbox_inches="tight",
        dpi=300,
    )


def get_reduced_ancestral_dict(ancestral: dict, src: str):
    # Calculate mut and missing
    mut = sum(
        ancestral.get(category, 0) for category in Constants.ANCESTRAL_CATEGORY["mut"]
    )
    missing = sum(
        ancestral.get(category, 0)
        for category in Constants.ANCESTRAL_CATEGORY["missing"]
    )

    # Calculate the ancestral_dict
    ancestral_dict = {
        key: value
        for key, value in ancestral.items()
        if key
        not in Constants.ANCESTRAL_CATEGORY["mut"]
        + Constants.ANCESTRAL_CATEGORY["missing"]
    }
    ancestral_dict["mut"] = mut
    ancestral_dict["missing"] = Constants.ANCESTRAL_NGENES[src] - sum(
        ancestral_dict.values()
    )  # 18430 IS CONSTANT

    return ancestral_dict


def make_scatter_for_mammals(path: str, ancestral: dict, source: str):
    fig, ax = plt.subplots()

    df = pd.read_csv(
        Constants.FileNames.MAMMALS,
        sep="\t",
        thousands=",",
    )
    df["mut"] = df[["Lost", "Uncertain loss"]].sum(axis=1)
    df["missing"] = df[["Missing", "Partially intact"]].sum(axis=1)
    df["taxonomy"] = [
        lineage.split("; ")[4]
        if lineage.split("; ")[4] in Constants.SUPERORDER
        else "Other"
        for lineage in df["Taxonomic Lineage"]
    ]
    df = df[["Species", "taxonomy", "Intact", "mut", "missing"]]

    ancestral_dict = get_reduced_ancestral_dict(ancestral, source)

    # # FAST FIX:
    # # if no genes were classified as "Intact", then number of keys in ancestral_dict will be two instead of three.
    # # As result, adding row to the dataframe will fail.
    # # Such story seems to happen in 100% of cases, if non ensembl gene ids are used.
    # # I tried to trace upstream, but stopped at qual_by_ancestral function from assembly_stats.py.
    # # Issue is a bit upstream of it
    #
    # ancestral_values = list(ancestral_dict.values())
    # if len(ancestral_values) == 2:
    #     ancestral_values = [0] + ancestral_values
    # df.loc[len(df.index), :] = ["User", "User"] + ancestral_values

    df.loc[len(df.index), :] = ["User", "User"] + list(ancestral_dict.values())
    df.iloc[:, 2:] = df.iloc[:, 2:] / Constants.ANCESTRAL_NGENES[source]

    for i, row in df.iterrows():
        if row.taxonomy != "User":
            ax.scatter(
                row.missing,
                row.mut,
                c=Constants.SUPERORDER_COLORS[row.taxonomy],
                alpha=0.4,
                label=row.taxonomy,
                marker="+",
            )
        else:
            ax.scatter(
                row.missing,
                row.mut,
                c=Constants.SUPERORDER_COLORS[row.taxonomy],
                alpha=1,
                label="Your assembly",
            )

    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    l = ax.legend(by_label.values(), by_label.keys(), frameon=False, fontsize="small")
    for text in l.get_texts():
        if text.get_text() == "Your assembly":
            text.set_color("red")
    ax.set_title("Mammalia: 510 species", fontsize=11, ha="center")
    ax.set_xlabel("%ancestral genes with\nmissing sequence")
    ax.set_ylabel("%ancestral genes with\ninactivating mutations")

    plt.savefig(
        os.path.join(
            path, Constants.DirNames.FIGURES, Constants.FigNames.ANCESTRAL_SCATTER
        ),
        format=Constants.FigNames.PLOTTING_FORMAT,
        bbox_inches="tight",
        dpi=300,
    )


def make_ancestral_barplot(path: str, ancestral: dict, source):
    fig, ax = plt.subplots(figsize=(7, 5))

    ancestral_dict = get_reduced_ancestral_dict(ancestral, source)
    ancestral_dict = {k: v / 18430 * 100 for k, v in ancestral_dict.items()}
    for x, y in ancestral_dict.items():
        ax.bar(x, y)

    bar_setup(ax)

    for i, col in enumerate(["Intact", "Inactivated", "Missing"]):
        ax.text(i, -6, col, ha="center", va="center", fontsize=10)

    for i, val in enumerate(ancestral_dict.values()):
        ax.text(i, val + 3, f"{val:.2f}%", ha="center", va="center", fontsize=10)

    ax.set_ylabel("%ancestral classes in query")
    ax.set_title("Ancestral classes in annotation", fontsize=11, ha="center")

    plt.savefig(
        os.path.join(
            path, Constants.DirNames.FIGURES, Constants.FigNames.ANCESTRAL_BARPLOT
        ),
        format=Constants.FigNames.PLOTTING_FORMAT,
        bbox_inches="tight",
        dpi=300,
    )


def make_lengths_histogram(path: str, lengths: pd.Series):
    fig, ax = plt.subplots(figsize=(7, 5))
    mean = np.mean(lengths)
    sd = np.std(lengths)
    median = np.median(lengths)

    n, bins, patches = ax.hist(
        lengths, bins=Constants.HIST_NBINS, density=True, alpha=0.5
    )
    n.sort()
    ax.text(bins[-20], n[-5], f"Mean: {mean:.6}\nMedian:{median:.6}\nStd: {sd:.6}")
    ax.set_xlabel("Gene lengths [bp]")
    ax.set_ylabel("Normalized amplitude")
    ax.set_title("Gene lengths distribution", fontsize=11, ha="center")

    plt.savefig(
        os.path.join(path, Constants.DirNames.FIGURES, Constants.FigNames.LENGTHS_PLOT),
        format=Constants.FigNames.PLOTTING_FORMAT,
        bbox_inches="tight",
        dpi=300,
    )


def make_scores_histogram(path: str, table: pd.DataFrame):
    fig, ax = plt.subplots(figsize=(7, 4))

    scores = table[table["pred"] > 0]["pred"]
    mean = np.mean(scores)
    sd = np.std(scores)
    median = np.median(scores)

    n, bins, patches = ax.hist(scores, bins=50, density=False, alpha=0.5)
    n.sort()
    ax.text(bins[5], n[-1] / 1.5, f"Mean: {mean:.6}\nMedian:{median:.6}\nStd: {sd:.6}")
    ax.set_xlabel("Prediction scores")
    ax.set_ylabel("Number of transcripts")
    ax.set_title("Orthology score distribution", fontsize=12, ha="center")

    plt.savefig(
        os.path.join(path, Constants.DirNames.FIGURES, Constants.FigNames.SCORE_PLOT),
        format=Constants.FigNames.PLOTTING_FORMAT,
        bbox_inches="tight",
        dpi=300,
    )


def make_busco_barplot(path: str, stats: list):
    ax = pd.DataFrame(stats).plot(kind="bar")
    ax.set_title("BUSCO completeness of your assembly")
    ax.set_ylabel("Completeness proportion (%)")

    for idx, pair in enumerate(stats):
        ax.text(idx - 0.15, 5, f"{pair[1]:.4}")

    ax.set_xticklabels([x[0] for x in stats], rotation=45)
    ax.legend().remove()
    ax.set_ylim(0, 100.5)

    plt.savefig(
        os.path.join(
            path, Constants.DirNames.FIGURES, Constants.FigNames.BUSCO_BARPLOT
        ),
        format=Constants.FigNames.PLOTTING_FORMAT,
        bbox_inches="tight",
        dpi=300,
    )


def postoga_plotter(
    path: str,
    table: pd.DataFrame,
    ancestral: dict,
    db: list,
    ngenes,
    lengths: pd.Series,
    src: str,
    busco_stats: list,
    species="human",
    db1=None,
):
    """The master plotter of postoga"""

    log = Log.connect(path, Constants.FileNames.LOG)
    log.record(f"postoga finished processing data! building report...")

    set_font()
    os.makedirs(os.path.join(path, Constants.DirNames.FIGURES), exist_ok=True)

    # fig = plt.figure(figsize=(7, 12))
    # plt.subplots_adjust(hspace=0.8)
    # gs = gridspec.GridSpec(4, 2, height_ratios=[0.09, 1, 1, 1])
    #
    # display metadata on plot
    # commit = shell(Constants.Commands.COMMIT)
    # branch = shell(Constants.Commands.BRANCH)

    # fig.text(0.67, 0.82, Constants.PLOTSTAMP.format(timestamp, __version__, branch, commit), fontsize=8, ha="center")
    #
    # fig.text(0.1, 0.79, "A", fontsize=15, ha="center")
    # fig.text(0.1, 0.525, "C", fontsize=15, ha="center")
    # fig.text(0.1, 0.28, "D", fontsize=15, ha="center")
    # fig.text(0.52, 0.28, "E", fontsize=15, ha="center")
    #
    # # display logo
    # logo_img = plt.imread(Constants.FileNames.LOGO_IMG) # "./supply/postoga_logo.png")
    # logo = OffsetImage(logo_img, zoom=0.2)
    # ax_logo = fig.add_subplot(gs[0, 0])
    # ax_logo.add_artist(AnnotationBbox(logo, (0.25, 0.75), frameon=False))
    # ax_logo.axis("off")

    # # init subplots
    # ax1 = fig.add_subplot(gs[1, 0])
    # ax2 = fig.add_subplot(gs[1, 1])
    # ax3 = fig.add_subplot(gs[2, :])
    # ax4 = fig.add_subplot(gs[3, 0])
    # ax5 = fig.add_subplot(gs[3, 1])
    #

    if db1:
        make_query_barplot_stats(path, [db, db1])
    else:
        make_query_barplot_stats(path, [db])

    make_scores_histogram(path, table)
    make_boxplot_annotated_genes(path, species, ngenes)
    make_ancestral_barplot(path, ancestral, src)
    make_scatter_for_mammals(path, ancestral, src)
    make_lengths_histogram(path, lengths)
    make_busco_barplot(path, busco_stats)
