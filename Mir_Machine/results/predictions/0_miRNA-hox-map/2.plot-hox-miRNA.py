#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import argparse
import re
import numpy as np
#RUN as: python3 3.plot.py --hox hox_clusters.annotated.csv --mir Canonical_HOX_miRNAs.tsv  --mirnas Mir-10.PRE Iab-4.PRE Mir-279.PRE Mir-2.PRE
############################################################
# ARGUMENTS
############################################################

parser = argparse.ArgumentParser()

parser.add_argument(
    "--hox",
    required=True,
    help="hox_clusters.annotated.csv"
)

parser.add_argument(
    "--mir",
    required=True,
    help="Canonical_HOX_miRNAs.tsv"
)

parser.add_argument(
    "--mirnas",
    nargs="+",
    default=[
        "Mir-10.PRE",
        "Iab-4.PRE"
    ]
)

args = parser.parse_args()

############################################################
# READ DATA
############################################################

hox = pd.read_csv(args.hox)

mirs = pd.read_csv(
    args.mir,
    sep="\t"
)

############################################################
# CLEAN TAXON NAMES
############################################################

def clean_taxon(x):

    x = str(x)

    x = re.sub(r"\.fasta$", "", x)
    x = re.sub(r"\.fa$", "", x)

    return x

hox["Taxon"] = hox["genome"].apply(
    clean_taxon
)

############################################################
# FILTER MIRNAS
############################################################

mirs = mirs[
    mirs["miRNA"].isin(
        set(args.mirnas)
    )
].copy()

mirs = mirs[
    mirs["Status"].isin(
        [
            "INSIDE_HOX",
            "SAME_SCAFFOLD"
        ]
    )
].copy()

############################################################
# TAXA
############################################################

taxa = sorted(
    set(hox["Taxon"])
    &
    set(mirs["Taxon"])
)

############################################################
# COUNT TRACKS
############################################################

n_tracks = 0

for taxon in taxa:

    n_tracks += len(
        hox[
            hox["Taxon"] == taxon
        ]["scaffold"].unique()
    )

fig_height = max(
    8,
    n_tracks * 0.5
)

fig, ax = plt.subplots(
    figsize=(16, fig_height)
)

############################################################
# MARKERS
############################################################

marker_dict = {

    "Mir-10.PRE": "o",

    "Iab-4.PRE": "s"

}

############################################################
# PLOT
############################################################

y = n_tracks + len(taxa)

for taxon in taxa:

    ax.text(

        -0.05,

        y,

        taxon,

        fontsize=12,

        fontweight="bold",

        ha="right",

        va="center",

        transform=ax.get_yaxis_transform()

    )

    tax_hox = hox[
        hox["Taxon"] == taxon
    ]

    scaffolds = sorted(
        tax_hox["scaffold"].unique()
    )

    for scaffold in scaffolds:

        y -= 1

        sdf = tax_hox[
            tax_hox["scaffold"] == scaffold
        ].copy()

        scaffold_mirs = mirs[
            (mirs["Taxon"] == taxon)
            &
            (mirs["Scaffold"] == scaffold)
        ].copy()

        ####################################################
        # COORDINATE RANGE
        ####################################################

        coords = []

        coords.extend(
            sdf["start"].tolist()
        )

        coords.extend(
            sdf["end"].tolist()
        )

        if len(scaffold_mirs) > 0:

            coords.extend(
                scaffold_mirs["Position"]
                .tolist()
            )

        xmin = min(coords)
        xmax = max(coords)

        span = xmax - xmin

        if span == 0:
            span = 1

        ####################################################
        # SCAFFOLD LINE
        ####################################################

        ax.plot(
            [0, 1],
            [y, y],
            lw=1.5,
            color="black"
        )

        ax.text(

            -0.01,

            y,

            scaffold,

            fontsize=7,

            ha="right",

            va="center"

        )

        ####################################################
        # HOX GENES
        ####################################################

        for _, gene in sdf.iterrows():

            start = (
                gene["start"] - xmin
            ) / span

            end = (
                gene["end"] - xmin
            ) / span

            width = max(
                0.005,
                end - start
            )

            rect = patches.Rectangle(

                (
                    start,
                    y - 0.15
                ),

                width,

                0.3

            )

            ax.add_patch(
                rect
            )

            ax.text(

                start + width/2,

                y + 0.22,

                gene["hoxgene"],

                fontsize=7,

                rotation=90,

                ha="center",

                va="bottom"

            )

        ####################################################
        # miRNAs
        ####################################################

        for _, mir in scaffold_mirs.iterrows():

            pos = (
                mir["Position"] - xmin
            ) / span

            marker = marker_dict.get(
                mir["miRNA"],
                "^"
            )

            if mir["Status"] == "INSIDE_HOX":

                yy = y

            else:

                yy = y + 0.25

            ax.scatter(

                pos,

                yy,

                marker=marker,

                s=80,

                zorder=10

            )

            ax.text(

                pos,

                yy + 0.15,

                mir["miRNA"]
                .replace(
                    ".PRE",
                    ""
                ),

                fontsize=6,

                rotation=45,

                ha="center"

            )

    y -= 1.2

############################################################
# LEGEND
############################################################

handles = []

for mirna in sorted(
    set(mirs["miRNA"])
):

    handles.append(

        plt.Line2D(

            [0],

            [0],

            marker=marker_dict.get(
                mirna,
                "^"
            ),

            linestyle="",

            markersize=8,

            label=mirna.replace(
                ".PRE",
                ""
            )

        )

    )

ax.legend(

    handles=handles,

    loc="upper left",

    bbox_to_anchor=(1.01, 1)

)

############################################################
# FORMATTING
############################################################

ax.set_xlim(
    -0.05,
    1.02
)

ax.set_ylim(
    0,
    n_tracks + len(taxa) + 2
)

ax.set_xticks([])

ax.set_yticks([])

ax.set_xlabel(
    "Relative position along scaffold"
)

ax.set_title(
    "Genomic positions of Hox genes and Hox-associated miRNAs"
)

plt.tight_layout()

############################################################
# OUTPUT
############################################################

plt.savefig(
    "PanelA_Hox_miRNA_map.pdf",
    bbox_inches="tight"
)

plt.savefig(
    "PanelA_Hox_miRNA_map.svg",
    bbox_inches="tight"
)

plt.close()

print("Written:")
print("  PanelA_Hox_miRNA_map.pdf")
print("  PanelA_Hox_miRNA_map.svg")
