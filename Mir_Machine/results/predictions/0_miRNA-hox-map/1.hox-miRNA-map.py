#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
import re

############################################################
# ARGUMENTS
#
#RUN as:  python3 1.hox-miRNA-map.py -m gffs_hox-miRNA --hox ./hox_clusters.annotated.csv
############################################################

parser = argparse.ArgumentParser()

parser.add_argument(
    "-m",
    "--mirna_dir",
    required=True,
    help="Directory containing *_miRNA.gff"
)

parser.add_argument(
    "--hox",
    required=True,
    help="HOX coordinate csv"
)

args = parser.parse_args()

############################################################
# READ HOX CSV
############################################################

hox = pd.read_csv(args.hox)

required = [
    "genome",
    "hoxgene",
    "scaffold",
    "start",
    "end"
]

for col in required:

    if col not in hox.columns:

        raise ValueError(
            f"Missing column '{col}' in HOX file"
        )

############################################################
# TAXON NAMES
############################################################

def get_taxon(genome):

    genome = str(genome)

    genome = re.sub(r"\.fasta$", "", genome)
    genome = re.sub(r"\.fa$", "", genome)

    return genome

hox["Taxon"] = hox["genome"].apply(
    get_taxon
)

############################################################
# BUILD HOX CLUSTERS
############################################################

clusters = {}

for taxon, tdf in hox.groupby("Taxon"):

    for scaffold, sdf in tdf.groupby("scaffold"):

        clusters[(taxon, scaffold)] = {

            "start":
                sdf["start"].min(),

            "end":
                sdf["end"].max(),

            "genes":
                list(
                    sdf["hoxgene"]
                )
        }

############################################################
# ORDERED HOX MAP
############################################################

hox["midpoint"] = (
    hox["start"] +
    hox["end"]
) / 2

hox_order = {}

for (taxon, scaffold), sdf in hox.groupby(
    ["Taxon", "scaffold"]
):

    hox_order[(taxon, scaffold)] = (
        sdf
        .sort_values("midpoint")
        .reset_index(drop=True)
    )

############################################################
# READ MIRMACHINE GFF
############################################################

def read_mir_gff(gff):

    rows = []

    with open(gff) as f:

        for line in f:

            if line.startswith("#"):
                continue

            cols = line.rstrip().split("\t")

            if len(cols) < 9:
                continue

            scaffold = cols[0]

            start = int(cols[3])
            end   = int(cols[4])

            attr = cols[8]

            m = re.search(
                r'gene_id=([^;]+)',
                attr
            )

            family = (
                m.group(1)
                if m
                else "Unknown"
            )

            rows.append([

                scaffold,

                start,

                end,

                family

            ])

    return pd.DataFrame(

        rows,

        columns=[

            "scaffold",

            "start",

            "end",

            "family"

        ]

    )

############################################################
# PROCESS TAXA
############################################################

results = []

mir_dir = Path(args.mirna_dir)

for mir_file in sorted(
    mir_dir.glob("*_miRNA.gff")
):

    taxon = mir_file.name.replace(
        "_miRNA.gff",
        ""
    )

    print(
        f"Processing {taxon}"
    )

    mirs = read_mir_gff(
        mir_file
    )

    for _, mir in mirs.iterrows():

        center = (
            mir["start"] +
            mir["end"]
        ) / 2

        scaffold = mir["scaffold"]

        context = "NO_HOX_LINKAGE"

        ####################################################
        # SAME SCAFFOLD AS HOX?
        ####################################################

        if (taxon, scaffold) in clusters:

            cluster = clusters[
                (taxon, scaffold)
            ]

            cluster_start = cluster["start"]
            cluster_end   = cluster["end"]

            ################################################
            # DISTANCE TO CLUSTER
            ################################################

            if (
                cluster_start
                <= center
                <= cluster_end
            ):

                distance = 0

                status = "INSIDE_HOX"

            elif center < cluster_start:

                distance = (
                    cluster_start
                    - center
                )

                status = (
                    "SAME_SCAFFOLD"
                )

            else:

                distance = (
                    center
                    - cluster_end
                )

                status = (
                    "SAME_SCAFFOLD"
                )

            ################################################
            # HOX CONTEXT
            ################################################

            genes = hox_order[
                (taxon, scaffold)
            ]

            found = False

            for _, gene in genes.iterrows():

                if (
                    gene["start"]
                    <= center
                    <= gene["end"]
                ):

                    context = (
                        f"inside_"
                        f"{gene['hoxgene']}"
                    )

                    found = True

                    break

            if not found:

                mids = genes[
                    "midpoint"
                ].values

                names = genes[
                    "hoxgene"
                ].values

                if center < mids[0]:

                    context = (
                        f"upstream_of_"
                        f"{names[0]}"
                    )

                elif center > mids[-1]:

                    context = (
                        f"downstream_of_"
                        f"{names[-1]}"
                    )

                else:

                    for i in range(
                        len(mids) - 1
                    ):

                        if (
                            mids[i]
                            <= center
                            <= mids[i+1]
                        ):

                            context = (
                                f"between_"
                                f"{names[i]}"
                                f"_and_"
                                f"{names[i+1]}"
                            )

                            break

        ####################################################
        # DIFFERENT SCAFFOLD
        ####################################################

        else:

            distance = np.nan

            status = (
                "DIFFERENT_SCAFFOLD"
            )

        ####################################################
        # SCORE
        ####################################################

        if status == "DIFFERENT_SCAFFOLD":

            score = 0

        else:

            score = (
                10 -
                np.log10(
                    distance + 1
                )
            )

        results.append([

            taxon,

            mir["family"],

            scaffold,

            int(center),

            status,

            distance,

            score,

            context

        ])

############################################################
# OUTPUT TABLE
############################################################

out = pd.DataFrame(

    results,

    columns=[

        "Taxon",

        "miRNA",

        "Scaffold",

        "Position",

        "Status",

        "Distance_bp",

        "HOX_score",

        "HOX_context"

    ]

)

out.to_csv(

    "miRNA_HOX_positions.tsv",

    sep="\t",

    index=False

)

############################################################
# PRESENCE ABSENCE
############################################################

presence = (
    out.assign(
        present=1
    )
    .pivot_table(
        index="Taxon",
        columns="miRNA",
        values="present",
        aggfunc="max",
        fill_value=0
    )
)

presence.to_csv(
    "miRNA_presence_absence.tsv",
    sep="\t"
)

############################################################
# CANONICAL HOX miRNAs
############################################################

hox_mirs = out[

    out["miRNA"].str.contains(

        "Mir-10|Iab-4|Mir-2",

        case=False,

        na=False

    )

].sort_values(

    ["Taxon", "HOX_score"],

    ascending=[True, False]

)

hox_mirs.to_csv(

    "Canonical_HOX_miRNAs.tsv",

    sep="\t",

    index=False

)

############################################################
# HEATMAP MATRIX
############################################################

matrix = out.pivot_table(

    index="Taxon",

    columns="miRNA",

    values="HOX_score",

    aggfunc="max",

    fill_value=0

)

matrix.to_csv(

    "All_miRNA_HOX_distance_matrix.tsv",

    sep="\t"

)

############################################################
# HEATMAP
############################################################

fig, ax = plt.subplots(

    figsize=(

        max(
            10,
            matrix.shape[1] * 0.35
        ),

        max(
            6,
            matrix.shape[0] * 0.4
        )

    )

)

im = ax.imshow(

    matrix,

    aspect="auto",

    cmap="Reds",

    vmin=0,

    vmax=10

)

ax.set_xticks(
    np.arange(
        matrix.shape[1]
    )
)

ax.set_xticklabels(

    matrix.columns,

    rotation=90,

    fontsize=7

)

ax.set_yticks(
    np.arange(
        matrix.shape[0]
    )
)

ax.set_yticklabels(

    matrix.index,

    fontsize=8

)

cbar = plt.colorbar(im)

cbar.set_label(
    "HOX association score (10 − log10(distance+1))"
)

plt.tight_layout()

plt.savefig(
    "miRNA_HOX_heatmap.pdf"
)

plt.savefig(
    "miRNA_HOX_heatmap.svg"
)

plt.close()

############################################################
# FINISHED
############################################################

print("\nOutput:")

print("  miRNA_HOX_positions.tsv")
print("  miRNA_presence_absence.tsv")
print("  Canonical_HOX_miRNAs.tsv")
print("  All_miRNA_HOX_distance_matrix.tsv")
print("  miRNA_HOX_heatmap.pdf")
print("  miRNA_HOX_heatmap.svg")
