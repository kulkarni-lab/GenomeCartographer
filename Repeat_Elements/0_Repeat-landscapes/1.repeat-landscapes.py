#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import argparse
from scipy.signal import savgol_filter

###########################################################
# ARGUMENTS
###########################################################

# RUN:
# python3 1.repeat-landscapes.py -i ./out_files/ -g gen-sizes.txt

parser = argparse.ArgumentParser()

parser.add_argument(
    "-i",
    "--indir",
    required=True,
    help="Directory containing RepeatMasker .out files"
)

parser.add_argument(
    "-g",
    "--genomes",
    required=True,
    help="Species genome size table"
)

args = parser.parse_args()

###########################################################
# GENOME SIZES
###########################################################

genome_sizes = pd.read_csv(
    args.genomes,
    sep=r"\s+"
)

genome_sizes = dict(
    zip(
        genome_sizes["species"],
        genome_sizes["genome_size"]
    )
)

###########################################################
# HELPER FUNCTIONS
###########################################################

def simplify_class(x):

    if x.startswith("LINE"):
        return "LINE"

    elif x.startswith("LTR"):
        return "LTR"

    elif x.startswith("SINE"):
        return "SINE"

    elif x.startswith("DNA"):
        return "DNA"

    elif x.startswith("RC"):
        return "Helitron"

    else:
        return "Unknown"


def smooth_curve(y, window=7, poly=3):

    y = np.asarray(y)

    if len(y) < window:
        return y

    if window % 2 == 0:
        window += 1

    return savgol_filter(
        y,
        window_length=window,
        polyorder=poly
    )


def weighted_quantile(values, weights, quantile):

    values = np.asarray(values)
    weights = np.asarray(weights)

    order = np.argsort(values)

    values = values[order]
    weights = weights[order]

    cumulative_weights = np.cumsum(weights)

    cutoff = quantile * weights.sum()

    idx = np.searchsorted(
        cumulative_weights,
        cutoff
    )

    return values[idx]


def weighted_median(values, weights):

    return weighted_quantile(
        values,
        weights,
        0.5
    )


def load_repeatmasker_out(outfile):

    rows = []

    with open(outfile) as f:

        for line in f:

            if (
                line.startswith("SW")
                or line.startswith("score")
                or line.strip() == ""
            ):
                continue

            parts = line.split()

            if len(parts) < 11:
                continue

            try:

                div = float(parts[1])

                start = int(parts[5])
                end = int(parts[6])

                length = end - start + 1

                rclass = parts[10]

            except Exception:
                continue

            rows.append([
                div,
                length,
                rclass
            ])

    return pd.DataFrame(
        rows,
        columns=[
            "div",
            "length",
            "class"
        ]
    )

###########################################################
# PROCESS FILES
###########################################################

indir = Path(args.indir)

all_species = {}

bins = np.arange(0, 41, 1)

div_summary = []

for outfile in sorted(indir.glob("*.out")):

    species = outfile.stem

    if species not in genome_sizes:

        print(
            f"Skipping {species}: no genome size"
        )
        continue

    genome_size = genome_sizes[species]

    print(f"\nProcessing {species}")

    df = load_repeatmasker_out(outfile)

    if len(df) == 0:

        print("No repeat entries found")

        continue

    df["class_simple"] = (
        df["class"]
        .astype(str)
        .apply(simplify_class)
    )

    #######################################################
    # DIVERGENCE STATISTICS
    #######################################################

    weighted_mean_div = np.average(
        df["div"],
        weights=df["length"]
    )

    weighted_sd_div = np.sqrt(
        np.average(
            (df["div"] - weighted_mean_div) ** 2,
            weights=df["length"]
        )
    )

    median_div = weighted_median(
        df["div"],
        df["length"]
    )

    q25_div = weighted_quantile(
        df["div"],
        df["length"],
        0.25
    )

    q75_div = weighted_quantile(
        df["div"],
        df["length"],
        0.75
    )

    iqr_div = q75_div - q25_div

    print(
        f"  Mean divergence   : {weighted_mean_div:.2f}%"
    )

    print(
        f"  Median divergence : {median_div:.2f}%"
    )

    print(
        f"  Q1 divergence     : {q25_div:.2f}%"
    )

    print(
        f"  Q3 divergence     : {q75_div:.2f}%"
    )

    print(
        f"  IQR divergence    : {iqr_div:.2f}"
    )

    print(
        f"  SD divergence     : {weighted_sd_div:.2f}"
    )

    div_summary.append([
        species,
        weighted_mean_div,
        median_div,
        q25_div,
        q75_div,
        iqr_div,
        weighted_sd_div
    ])

    all_species[species] = (
        df,
        genome_size
    )

###########################################################
# INDIVIDUAL LANDSCAPES
###########################################################

for species, (df, genome_size) in all_species.items():

    plt.figure(figsize=(8, 5))

    for te_class in [
        "DNA",
        "LINE",
        "LTR",
        "SINE",
        "Helitron",
        "Unknown"
    ]:

        sub = df[
            df["class_simple"] == te_class
        ]

        if len(sub) == 0:
            continue

        percents = []

        for i in range(len(bins) - 1):

            lo = bins[i]
            hi = bins[i + 1]

            bp = sub[
                (sub["div"] >= lo)
                &
                (sub["div"] < hi)
            ]["length"].sum()

            perc = (
                bp /
                genome_size
            ) * 100

            percents.append(perc)

        x = bins[:-1] + 0.5

        y = smooth_curve(
            percents,
            window=7,
            poly=3
        )

        plt.plot(
            x,
            y,
            linewidth=2,
            label=te_class
        )

    plt.xlabel(
        "Percent divergence from consensus"
    )

    plt.ylabel(
        "% genome occupied"
    )

    plt.title(species)

    plt.legend(
        fontsize=8
    )

    plt.tight_layout()

    plt.savefig(
        f"{species}_landscape_percentGenome.pdf"
    )

    plt.savefig(
        f"{species}_landscape_percentGenome.svg"
    )

    plt.close()

###########################################################
# COMBINED FIGURE
###########################################################

plt.figure(figsize=(10, 6))

for species, (df, genome_size) in all_species.items():

    percents = []

    for i in range(len(bins) - 1):

        lo = bins[i]
        hi = bins[i + 1]

        bp = df[
            (df["div"] >= lo)
            &
            (df["div"] < hi)
        ]["length"].sum()

        perc = (
            bp /
            genome_size
        ) * 100

        percents.append(perc)

    peak_idx = np.argmax(percents)

    peak_divergence = (
        bins[:-1] + 0.5
    )[peak_idx]

    peak_occupancy = percents[peak_idx]

    print(
        f"  Peak divergence   : "
        f"{peak_divergence:.1f}%"
    )

    print(
        f"  Peak occupancy    : "
        f"{peak_occupancy:.4f}% genome"
    )

    x = bins[:-1] + 0.5

    y = smooth_curve(
        percents,
        window=7,
        poly=3
    )

    plt.plot(
        x,
        y,
        linewidth=2,
        label=species
    )

plt.xlabel(
    "Percent divergence from consensus"
)

plt.ylabel(
    "% genome occupied"
)

plt.title(
    "TE landscapes across species"
)

plt.legend(
    bbox_to_anchor=(1.02, 1),
    loc="upper left",
    fontsize=8
)

plt.tight_layout()

plt.savefig(
    "All_species_TE_landscape_percentGenome.pdf"
)

plt.savefig(
    "All_species_TE_landscape_percentGenome.svg"
)

plt.close()

###########################################################
# SAVE DIVERGENCE SUMMARY
###########################################################

summary_df = pd.DataFrame(
    div_summary,
    columns=[
        "Species",
        "Mean_divergence",
        "Median_divergence",
        "Q1_divergence",
        "Q3_divergence",
        "IQR_divergence",
        "SD_divergence"
    ]
)

summary_df.to_csv(
    "TE_divergence_summary.tsv",
    sep="\t",
    index=False
)

print(
    "\nSaved TE_divergence_summary.tsv"
)

