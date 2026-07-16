#!/usr/bin/env python3
"""
Genome architecture analysis
"""

import os
import re
import sys
import glob
import time
from collections import defaultdict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ranksums
from statsmodels.stats.multitest import multipletests

plt.rcParams["svg.fonttype"] = "none"

VIOLIN_COLOR = "#8d9091"
BOX_COLORS = [
    "#d95f5f", "#b39b2e", "#22b14c", "#12b3b6", "#7aa6e8",
    "#cc79a7", "#9b59b6", "#f39c12", "#95a5a6", "#34495e"
]


def parse_attrs(s):
    d = {}

    for k, v in re.findall(r'(\S+)\s+"([^"]+)"', s):
        d[k] = v

    if not d:
        for x in s.split(";"):
            if "=" in x:
                k, v = x.split("=", 1)
                d[k.strip()] = v.strip()

    return d


def load_order(f):
    with open(f, "r") as fh:
        return [x.strip() for x in fh if x.strip()]


def add_sig(ax, df, col, order, ref):

    ps = []
    comps = []

    rv = df[df.species == ref][col]

    for s in order:
        if s == ref:
            continue

        vv = df[df.species == s][col]

        if len(rv) < 2 or len(vv) < 2:
            continue

        _, p = ranksums(rv, vv)

        ps.append(p)
        comps.append(s)

    if not ps:
        return

    adj = multipletests(ps, method="bonferroni")[1]

    y = df[col].max() * 1.3
    step = (
        df[col].max() * 0.15
        if ax.get_yscale() != "log"
        else df[col].max() * 0.5
    )

    for s, p in zip(comps, adj):

        if p >= 0.05:
            continue

        x1 = order.index(ref) + 1
        x2 = order.index(s) + 1

        ax.plot(
            [x1, x1, x2, x2],
            [y, y * 1.05, y * 1.05, y],
            c="black",
            lw=1
        )

        lab = (
            "****" if p < 1e-4
            else "***" if p < 1e-3
            else "**" if p < 1e-2
            else "*"
        )

        ax.text((x1 + x2) / 2, y * 1.08, lab, ha="center")

        y += step


def process_file(path):

    species = os.path.basename(path).rsplit(".", 1)[0]

    print(f"    Reading {species} ...", flush=True)

    tx_exons = defaultdict(list)
    tx_cds = defaultdict(list)
    tx_gene = {}

    with open(path) as fh:

        for line in fh:

            if line.startswith("#"):
                continue

            f = line.rstrip().split("\t")

            if len(f) != 9:
                continue

            feature = f[2]
            start = int(f[3])
            end = int(f[4])

            a = parse_attrs(f[8])

            tid = (
                a.get("transcript_id")
                or a.get("Parent")
                or a.get("ID")
            )

            gid = (
                a.get("gene_id")
                or a.get("gene")
                or a.get("Parent")
            )

            if feature in ("mRNA", "transcript", "gene"):

                tid = a.get("transcript_id") or a.get("ID")
                gid = a.get("gene_id") or a.get("Parent") or tid

                if tid:
                    tx_gene[tid] = gid

            elif feature == "exon" and tid:
                tx_exons[tid].append((start, end))

            elif feature == "CDS" and tid:
                tx_cds[tid].append((start, end))

    gene_best = {}

    for tx in set(list(tx_exons.keys()) + list(tx_cds.keys())):

        feats = sorted(tx_exons[tx] if tx_exons[tx] else tx_cds[tx])

        if not feats:
            continue

        total = sum(e - s + 1 for s, e in feats)

        gid = tx_gene.get(tx, tx)

        if gid not in gene_best or total > gene_best[gid][0]:
            gene_best[gid] = (total, tx, feats)

    introns = []
    exons = []
    compact = []
    genes = []

    for gid, (tot, tx, feats) in gene_best.items():

        span = (
            max(e for s, e in feats)
            - min(s for s, e in feats)
            + 1
        )

        comp = tot / span if span > 0 else np.nan

        compact.append([species, gid, comp])

        for s, e in feats:
            exons.append([species, gid, e - s + 1])

        nintr = 0

        for i in range(len(feats) - 1):

            intr = feats[i + 1][0] - feats[i][1] - 1

            if intr > 0:
                nintr += 1
                introns.append([species, gid, intr])

        genes.append(
            [species, gid, span, len(feats), nintr, tot, comp]
        )

    print(
        f"    {species}: "
        f"{len(gene_best):,} genes retained, "
        f"{len(exons):,} exons, "
        f"{len(introns):,} introns",
        flush=True
    )

    return introns, exons, compact, genes


if __name__ == "__main__":

    if len(sys.argv) < 4:
        print(
            "Usage: python genome_architecture.py "
            "<gtf_dir> <prefix> <order.txt>"
        )
        sys.exit(1)

    gtfdir = sys.argv[1]
    prefix = sys.argv[2]
    order_file = sys.argv[3]

    order = load_order(order_file)
    ref = order[0]

    ai = []
    ae = []
    ac = []
    ag = []

    files = (
        glob.glob(os.path.join(gtfdir, "*.gtf"))
        + glob.glob(os.path.join(gtfdir, "*.gff"))
        + glob.glob(os.path.join(gtfdir, "*.gff3"))
    )

    if not files:
        print(f"No annotation files found in {gtfdir}")
        sys.exit(1)

    print(
        f"Found {len(files)} annotation files.\n",
        flush=True
    )

    total_files = len(files)

    for idx, f in enumerate(files, start=1):

        species = os.path.basename(f).rsplit(".", 1)[0]

        print(
            f"[{idx}/{total_files}] STARTING: {species}",
            flush=True
        )

        start_time = time.time()

        try:

            i, e, c, g = process_file(f)

            ai.extend(i)
            ae.extend(e)
            ac.extend(c)
            ag.extend(g)

            elapsed = time.time() - start_time

            print(
                f"[{idx}/{total_files}] DONE: {species} | "
                f"genes={len(g):,} "
                f"exons={len(e):,} "
                f"introns={len(i):,} "
                f"time={elapsed:.1f}s",
                flush=True
            )

        except Exception as exc:

            elapsed = time.time() - start_time

            print(
                f"[{idx}/{total_files}] ERROR: {species} "
                f"after {elapsed:.1f}s",
                flush=True
            )

            print(f"Reason: {exc}", flush=True)

            raise

    print("\nAll taxa processed.\n", flush=True)

    intr = pd.DataFrame(
        ai,
        columns=["species", "gene", "intron_length"]
    )

    ex = pd.DataFrame(
        ae,
        columns=["species", "gene", "exon_length"]
    )

    comp = pd.DataFrame(
        ac,
        columns=["species", "gene", "compactness"]
    )

    genes = pd.DataFrame(
        ag,
        columns=[
            "species",
            "gene",
            "gene_span",
            "num_exons",
            "num_introns",
            "total_exon_length",
            "compactness"
        ]
    )

    print("Writing output tables...", flush=True)

    intr.to_csv(
        f"{prefix}_introns.tsv",
        sep="\t",
        index=False
    )

    ex.to_csv(
        f"{prefix}_exons.tsv",
        sep="\t",
        index=False
    )

    comp.to_csv(
        f"{prefix}_compactness.tsv",
        sep="\t",
        index=False
    )

    genes.to_csv(
        f"{prefix}_gene_summary.tsv",
        sep="\t",
        index=False
    )

    print("Output tables written.", flush=True)

print("Output tables written.", flush=True)

# ------------------------------------------------------------
# SPECIES SUMMARY
# ------------------------------------------------------------

summary_rows = []

for species in order:

    if species not in genes["species"].unique():
        continue

    introns_species = intr[
        intr["species"] == species
    ]["intron_length"]

    exons_species = ex[
        ex["species"] == species
    ]["exon_length"]

    compact_species = comp[
        comp["species"] == species
    ]["compactness"]

    genes_species = genes[
        genes["species"] == species
    ]

    summary_rows.append({

        "species": species,

        "median_intron_length":
            introns_species.median(),

        "mean_intron_length":
            introns_species.mean(),

        "median_exon_length":
            exons_species.median(),

        "mean_exon_length":
            exons_species.mean(),

        "median_compactness":
            compact_species.median(),

        "mean_compactness":
            compact_species.mean(),

        "mean_introns_per_gene":
            genes_species["num_introns"].mean(),

        "median_gene_span":
            genes_species["gene_span"].median(),

        "num_genes":
            len(genes_species)
    })

summary_df = pd.DataFrame(summary_rows)

summary_df.to_csv(
    f"{prefix}_species_summary.tsv",
    sep="\t",
    index=False
)

print("Species summary written.", flush=True)

# ------------------------------------------------------------
# INTRON VIOLIN PLOT
# ------------------------------------------------------------

species_present = [
    s for s in order
    if s in intr["species"].unique()
]

if len(species_present) > 0:

    fig, ax = plt.subplots(figsize=(12, 6))

    violin_data = [
        intr[intr["species"] == s]["intron_length"]
        for s in species_present
    ]

    vp = ax.violinplot(
        violin_data,
        showmeans=True,
        showextrema=False
    )

    for body in vp["bodies"]:
        body.set_facecolor(VIOLIN_COLOR)
        body.set_edgecolor("black")
        body.set_alpha(0.85)

    ax.set_xticks(
        range(1, len(species_present) + 1)
    )

    ax.set_xticklabels(
        species_present,
        rotation=45,
        ha="right"
    )

    ax.set_yscale("log")

    ax.set_ylabel("Intron length")
    ax.set_title("Intron Size Distributions")

    add_sig(
        ax,
        intr,
        "intron_length",
        species_present,
        ref
    )

    plt.tight_layout()

    plt.savefig(
        f"{prefix}_intron_violin.svg",
        format="svg"
    )

    plt.close()

    print("Intron violin plot written.", flush=True)

# ------------------------------------------------------------
# EXON BOXPLOT
# ------------------------------------------------------------

species_present = [
    s for s in order
    if s in ex["species"].unique()
]

if len(species_present) > 0:

    fig, ax = plt.subplots(figsize=(12, 6))

    data = [
        ex[
            ex["species"] == s
        ]["exon_length"]
        for s in species_present
    ]

    bp = ax.boxplot(
        data,
        showfliers=False,
        patch_artist=True
    )

    for i, patch in enumerate(bp["boxes"]):

        patch.set_facecolor(
            BOX_COLORS[
                i % len(BOX_COLORS)
            ]
        )

        patch.set_edgecolor("black")

    ax.set_xticks(
        range(1, len(species_present) + 1)
    )

    ax.set_xticklabels(
        species_present,
        rotation=45,
        ha="right"
    )

    ax.set_yscale("log")

    ax.set_ylabel("Exon length")
    ax.set_title("Exon Length Distribution")

    add_sig(
        ax,
        ex,
        "exon_length",
        species_present,
        ref
    )

    plt.tight_layout()

    plt.savefig(
        f"{prefix}_exon_boxplot.svg",
        format="svg"
    )

    plt.close()

    print("Exon boxplot written.", flush=True)

# ------------------------------------------------------------
# COMPACTNESS BOXPLOT
# ------------------------------------------------------------

species_present = [
    s for s in order
    if s in comp["species"].unique()
]

if len(species_present) > 0:

    fig, ax = plt.subplots(figsize=(12, 6))

    data = [
        comp[
            comp["species"] == s
        ]["compactness"]
        for s in species_present
    ]

    bp = ax.boxplot(
        data,
        showfliers=False,
        patch_artist=True
    )

    for i, patch in enumerate(bp["boxes"]):

        patch.set_facecolor(
            BOX_COLORS[
                i % len(BOX_COLORS)
            ]
        )

        patch.set_edgecolor("black")

    ax.set_xticks(
        range(1, len(species_present) + 1)
    )

    ax.set_xticklabels(
        species_present,
        rotation=45,
        ha="right"
    )

    ax.set_ylabel("Compactness")
    ax.set_ylim(0, 1.05)

    ax.set_title("Gene Compactness")

    add_sig(
        ax,
        comp,
        "compactness",
        species_present,
        ref
    )

    plt.tight_layout()

    plt.savefig(
        f"{prefix}_compactness_boxplot.svg",
        format="svg"
    )

    plt.close()

    print("Compactness boxplot written.", flush=True)

print()
print("Done.")
print(f"Outputs written with prefix: {prefix}")
print()
