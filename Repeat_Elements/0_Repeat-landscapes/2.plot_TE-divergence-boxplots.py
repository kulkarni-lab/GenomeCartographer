#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse

###########################################################
# ARGUMENTS
###########################################################

parser = argparse.ArgumentParser()

parser.add_argument(
    "-i",
    "--input",
    required=True,
    help="TE_divergence_summary.tsv"
)

args = parser.parse_args()

###########################################################
# LOAD DATA
###########################################################

df = pd.read_csv(
    args.input,
    sep="\t"
)

###########################################################
# SORT BY MEDIAN DIVERGENCE
###########################################################

df = df.sort_values(
    "Median_divergence"
).reset_index(drop=True)

species = df["Species"]

###########################################################
# FIGURE 1
# MEDIAN ± IQR (RECOMMENDED)
###########################################################

fig, ax = plt.subplots(
    figsize=(10, 6)
)

x = np.arange(len(df))

median = df["Median_divergence"]
q1 = df["Q1_divergence"]
q3 = df["Q3_divergence"]

lower = median - q1
upper = q3 - median

ax.errorbar(
    x,
    median,
    yerr=[lower, upper],
    fmt="o",
    capsize=5,
    linewidth=2
)

ax.set_xticks(x)

ax.set_xticklabels(
    species,
    rotation=60,
    ha="right"
)

ax.set_ylabel(
    "Divergence (%)"
)

ax.set_title(
    "Median TE divergence ± IQR"
)

plt.tight_layout()

plt.savefig(
    "Median_IQR_divergence.pdf"
)

plt.savefig(
    "Median_IQR_divergence.svg"
)

plt.close()

###########################################################
# FIGURE 2
# BOX-STYLE SUMMARY
###########################################################

fig, ax = plt.subplots(
    figsize=(11, 6)
)

for i, row in df.iterrows():

    q1 = row["Q1_divergence"]
    q3 = row["Q3_divergence"]
    median = row["Median_divergence"]
    sd = row["SD_divergence"]

    ax.add_patch(
        plt.Rectangle(
            (i - 0.3, q1),
            0.6,
            q3 - q1,
            fill=False
        )
    )

    ax.plot(
        [i - 0.3, i + 0.3],
        [median, median]
    )

    ax.plot(
        [i, i],
        [
            median - sd,
            median + sd
        ]
    )

ax.set_xticks(
    np.arange(len(df))
)

ax.set_xticklabels(
    species,
    rotation=60,
    ha="right"
)

ax.set_ylabel(
    "Divergence (%)"
)

ax.set_title(
    "TE divergence distribution summary"
)

plt.tight_layout()

plt.savefig(
    "TE_divergence_box_summary.pdf"
)

plt.savefig(
    "TE_divergence_box_summary.svg"
)

plt.close()

###########################################################
# FIGURE 3
# BUBBLE PLOT
###########################################################

fig, ax = plt.subplots(
    figsize=(8, 6)
)

bubble_sizes = (
    df["IQR_divergence"] * 120
)

scatter = ax.scatter(
    df["Median_divergence"],
    df["Mean_divergence"],
    s=bubble_sizes,
    alpha=0.7
)

for _, row in df.iterrows():

    ax.text(
        row["Median_divergence"],
        row["Mean_divergence"],
        row["Species"],
        fontsize=8
    )

ax.set_xlabel(
    "Median divergence (%)"
)

ax.set_ylabel(
    "Mean divergence (%)"
)

ax.set_title(
    "TE age structure across species"
)

plt.tight_layout()

plt.savefig(
    "TE_divergence_bubbleplot.pdf"
)

plt.savefig(
    "TE_divergence_bubbleplot.svg"
)

plt.close()

###########################################################
# FIGURE 4
# HEATMAP OF SUMMARY STATISTICS
###########################################################

stats = df[
    [
        "Mean_divergence",
        "Median_divergence",
        "Q1_divergence",
        "Q3_divergence",
        "IQR_divergence",
        "SD_divergence"
    ]
]

fig, ax = plt.subplots(
    figsize=(8, 8)
)

im = ax.imshow(
    stats,
    aspect="auto"
)

ax.set_xticks(
    np.arange(stats.shape[1])
)

ax.set_xticklabels(
    stats.columns,
    rotation=45,
    ha="right"
)

ax.set_yticks(
    np.arange(len(df))
)

ax.set_yticklabels(
    species
)

plt.colorbar(
    im,
    ax=ax,
    label="Divergence (%)"
)

ax.set_title(
    "Summary statistics heatmap"
)

plt.tight_layout()

plt.savefig(
    "TE_divergence_heatmap.pdf"
)

plt.savefig(
    "TE_divergence_heatmap.svg"
)

plt.close()

print(
    "Generated:"
)
print(
    "  Median_IQR_divergence.pdf/svg"
)
print(
    "  TE_divergence_box_summary.pdf/svg"
)
print(
    "  TE_divergence_bubbleplot.pdf/svg"
)
print(
    "  TE_divergence_heatmap.pdf/svg"
)
