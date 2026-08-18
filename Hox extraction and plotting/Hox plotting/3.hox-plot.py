#!/usr/bin/env python3

"""
plot_hox_clusters_from_chrom.py

RUN as: python3 3.hox-plot.py -c hox_genes.csv  -genomedir ../genomes/ -chromdir chroms/ -o HoxPlot_Indiac3

INPUT CSV FORMAT:
genome_name,chrom_file,sequence_name,hoxgene

Example:
Sarc.fna,Sarc.chrom,g1.t1,lab
Sarc.fna,Sarc.chrom,g55.t1,pb
Sarc.fna,Sarc.chrom,g212.t1,abdA

CHROM FILE FORMAT:
gene_id scaffold strand start end

Example:
g58.t1  CM075101.1_Phytoseiulus_persimilis_chromosome_1  +  4271    35893
g100.t1 NW_02323232.1_some_description                   -  132167  135632

OUTPUT:
    output/
        hox_clusters.annotated.csv
        hox_clusters.png
        hox_clusters.pdf

USAGE:
python plot_hox_clusters_from_chrom.py \
    -c hox_genes.csv \
    -genomedir ./genomes \
    -chromdir ./chroms \
    -o HoxPlot
"""

import argparse
import csv
import os
import re
import sys

import pandas as pd
import matplotlib.pyplot as plt

from matplotlib.patches import FancyArrowPatch


###############################################################################
# Colors
###############################################################################

HOX_COLORS = {

    "lab": "#d73027",
    "pb": "#fc8d59",
    "dfd": "#fee090",
    "scr": "#91cf60",
    "antp": "#1a9850",
    "ubx": "#4575b4",
    "abda": "#313695",
    "abdb": "#762a83",

    "post1": "#c51b7d",
    "post2": "#f0027f"
}


###############################################################################
# Normalize names
###############################################################################

def normalize_hox(h):

    h = h.lower()

    h = h.replace("-", "")
    h = h.replace("_", "")

    return h


###############################################################################
# Clean scaffold names
###############################################################################

def clean_scaffold(raw_scaffold):

    """
    Examples:
        CM075101.1_long_text      -> CM075101.1
        NW_02323232.1_extra       -> NW_02323232.1
        NC_000001.11_human_chr1   -> NC_000001.11
    """

    match = re.match(
        r'^([A-Z]{1,4}(?:_[0-9]+)?\.[0-9]+)',
        raw_scaffold
    )

    if match:
        return match.group(1)

    return raw_scaffold


###############################################################################
# Parse chrom file
###############################################################################

def parse_chrom_file(chrom_file):

    genes = {}

    print(f"[INFO] Parsing chrom file:")
    print(f"  {chrom_file}")

    with open(chrom_file) as f:

        for line in f:

            if line.startswith("#"):
                continue

            if line.strip() == "":
                continue

            fields = line.strip().split()

            if len(fields) < 5:
                continue

            gene_id = fields[0]

            ###################################################################
            # Clean scaffold name
            ###################################################################

            scaffold = clean_scaffold(fields[1])

            strand = fields[2]

            start = int(fields[3])
            end = int(fields[4])

            genes[gene_id] = {

                "scaffold": scaffold,
                "strand": strand,
                "start": start,
                "end": end

            }

    return genes


###############################################################################
# Plot clusters
###############################################################################

def plot_clusters(df, outfile_png, outfile_pdf):

    genomes = list(df["genome"].unique())

    fig_height = max(
        6,
        len(genomes) * 2.5
    )

    fig, ax = plt.subplots(
        figsize=(18, fig_height)
    )

    y_positions = {}

    for i, genome in enumerate(genomes[::-1]):

        y_positions[genome] = i

    ###########################################################################
    # Fixed display width for ALL taxa
    ###########################################################################

    DISPLAY_WIDTH = 1000000

    ###########################################################################
    # Draw genomes
    ###########################################################################

    for genome in genomes:

        subset = df[
            df["genome"] == genome
        ].copy()

        subset = subset.sort_values(
            ["scaffold", "start"]
        )

        y = y_positions[genome]

        #######################################################################
        # Genome label
        #######################################################################

        ax.text(
            -0.01,
            y,
            genome,
            transform=ax.get_yaxis_transform(),
            ha="right",
            va="center",
            fontsize=13,
            fontweight="bold"
        )

        #######################################################################
        # Compute scaffold lengths
        #######################################################################

        scaffold_offsets = {}

        offset = 0

        SCAFFOLD_SPACING = 150000

        total_length = 0

        scaffold_lengths = {}

        for scaffold in subset["scaffold"].unique():

            scaf_subset = subset[
                subset["scaffold"] == scaffold
            ]

            scaf_min = scaf_subset["start"].min()
            scaf_max = scaf_subset["end"].max()

            length = scaf_max - scaf_min

            scaffold_lengths[scaffold] = length

            total_length += length

        #######################################################################
        # Scaling factor
        #######################################################################

        scale_factor = DISPLAY_WIDTH / max(total_length, 1)

        #######################################################################
        # Assign scaffold offsets
        #######################################################################

        for scaffold in subset["scaffold"].unique():

            scaffold_offsets[scaffold] = offset

            scaled_len = (
                scaffold_lengths[scaffold]
                * scale_factor
            )

            offset += scaled_len + SCAFFOLD_SPACING

        #######################################################################
        # Draw scaffold baselines
        #######################################################################

        for scaffold in subset["scaffold"].unique():

            scaf = subset[
                subset["scaffold"] == scaffold
            ]

            xmin = scaffold_offsets[scaffold]

            xmax = (
                scaffold_offsets[scaffold]
                + scaffold_lengths[scaffold]
                * scale_factor
            )

            ax.plot(
                [xmin, xmax],
                [y, y],
                color="black",
                linewidth=1.5,
                alpha=0.6
            )

            ###################################################################
            # Scaffold labels
            ###################################################################

            ax.text(
                (xmin + xmax) / 2,
                y - 0.22,
                scaffold,
                ha="center",
                fontsize=7,
                alpha=0.8
            )

        #######################################################################
        # Draw genes
        #######################################################################

        MIN_GAP = 30000

        last_end = -999999999

        for _, row in subset.iterrows():

            scaf_subset = subset[
                subset["scaffold"] == row["scaffold"]
            ]

            scaf_min = scaf_subset["start"].min()

            start = (
                (
                    row["start"] - scaf_min
                )
                * scale_factor
                + scaffold_offsets[row["scaffold"]]
            )

            end = (
                (
                    row["end"] - scaf_min
                )
                * scale_factor
                + scaffold_offsets[row["scaffold"]]
            )

            ###################################################################
            # Avoid overlap
            ###################################################################

            if start - last_end < MIN_GAP:

                shift = MIN_GAP - (start - last_end)

                start += shift
                end += shift

            last_end = end

            gene = row["hoxgene"]

            color = HOX_COLORS.get(
                normalize_hox(gene),
                "gray"
            )

            ###################################################################
            # Reverse strand
            ###################################################################

            if row["strand"] == "-":

                arrow = FancyArrowPatch(

                    (end, y),
                    (start, y),

                    arrowstyle='Simple,tail_width=6,head_width=40,head_length=50',

                    linewidth=5,
                    color=color

                )

            ###################################################################
            # Forward strand
            ###################################################################

            else:

                arrow = FancyArrowPatch(

                    (start, y),
                    (end, y),

                    arrowstyle='Simple,tail_width=6,head_width=40,head_length=50',

                    linewidth=5,
                    color=color

                )

            ax.add_patch(arrow)

            ###################################################################
            # Gene labels
            ###################################################################

            ax.text(
                (start + end) / 2,
                y + 0.22,
                gene,
                ha="center",
                va="bottom",
                fontsize=10,
                fontweight="bold"
            )

    ###########################################################################
    # Remove axes
    ###########################################################################

    ax.set_yticks([])
    ax.set_xticks([])

    ax.set_xlabel("")

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)

    plt.tight_layout()

    plt.savefig(
        outfile_png,
        dpi=300,
        bbox_inches="tight"
    )

    plt.savefig(
        outfile_pdf,
        bbox_inches="tight"
    )

    print("\n[INFO] Figure written:")
    print(f"  - {outfile_png}")
    print(f"  - {outfile_pdf}")


###############################################################################
# Main
###############################################################################

def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-c",
        "--csv",
        required=True,
        help="Input CSV"
    )

    parser.add_argument(
        "-genomedir",
        "--genomedir",
        required=True,
        help="Genome FASTA directory"
    )

    parser.add_argument(
        "-chromdir",
        "--chromdir",
        required=True,
        help="Chrom file directory"
    )

    parser.add_argument(
        "-o",
        "--outdir",
        required=True,
        help="Output directory"
    )

    args = parser.parse_args()

    os.makedirs(
        args.outdir,
        exist_ok=True
    )

    ###########################################################################
    # Read CSV
    ###########################################################################

    rows = []

    with open(args.csv) as f:

        reader = csv.DictReader(f)

        for row in reader:
            rows.append(row)

    ###########################################################################
    # Validate files
    ###########################################################################

    for row in rows:

        genome_path = os.path.join(
            args.genomedir,
            row["genome_name"]
        )

        chrom_path = os.path.join(
            args.chromdir,
            row["chrom_file"]
        )

        if not os.path.isfile(genome_path):

            print(f"[ERROR] Missing genome:")
            print(f"  {genome_path}")

            sys.exit(1)

        if not os.path.isfile(chrom_path):

            print(f"[ERROR] Missing chrom file:")
            print(f"  {chrom_path}")

            sys.exit(1)

    ###########################################################################
    # Parse chrom files
    ###########################################################################

    chrom_cache = {}

    for row in rows:

        chrom_name = row["chrom_file"]

        if chrom_name not in chrom_cache:

            chrom_file = os.path.join(
                args.chromdir,
                chrom_name
            )

            chrom_cache[chrom_name] = parse_chrom_file(
                chrom_file
            )

    ###########################################################################
    # Annotate Hox genes
    ###########################################################################

    annotated = []

    for row in rows:

        genome = row["genome_name"]
        chrom_name = row["chrom_file"]
        seqid = row["sequence_name"]
        hox = row["hoxgene"]

        chrom_data = chrom_cache[chrom_name]

        found = False

        for gid, info in chrom_data.items():

            if gid == seqid or seqid in gid:

                annotated.append({

                    "genome": genome,
                    "chrom_file": chrom_name,
                    "sequence": seqid,
                    "hoxgene": hox,

                    "scaffold": info["scaffold"],
                    "start": info["start"],
                    "end": info["end"],
                    "strand": info["strand"]

                })

                found = True
                break

        if not found:

            print(f"[WARNING] Could not find:")
            print(f"  {seqid}")

    ###########################################################################
    # Write annotated CSV
    ###########################################################################

    outcsv = os.path.join(
        args.outdir,
        "hox_clusters.annotated.csv"
    )

    with open(outcsv, "w", newline="") as out:

        writer = csv.DictWriter(

            out,

            fieldnames=[

                "genome",
                "chrom_file",
                "sequence",
                "hoxgene",

                "scaffold",
                "start",
                "end",
                "strand"

            ]
        )

        writer.writeheader()

        for row in annotated:
            writer.writerow(row)

    print("\n[INFO] Annotated CSV written:")
    print(f"  - {outcsv}")

    ###########################################################################
    # Plot
    ###########################################################################

    df = pd.DataFrame(annotated)

    png = os.path.join(
        args.outdir,
        "hox_clusters.png"
    )

    pdf = os.path.join(
        args.outdir,
        "hox_clusters.pdf"
    )

    plot_clusters(
        df,
        png,
        pdf
    )


###############################################################################
# Run
###############################################################################

if __name__ == "__main__":
    main()
