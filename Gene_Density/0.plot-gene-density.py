#!/usr/bin/env python3

"""
Ridgeline (joyplot) of gene density distributions across genomes.

Input:
    gene_density/
        species1.txt
        species2.txt
        species3.txt
        ...

Each file must contain four tab-delimited columns:

    scaffold    start    end    density

Example:
    NW_003803370.1    0       100000    12.0
    NW_003803370.1    100000  200000    6.0

Output:
    gene_density_ridgeline.png
    gene_density_ridgeline.pdf
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

# --------------------------------------------------
# Input directory
# --------------------------------------------------

INPUT_DIR = "gene_density"

if not os.path.isdir(INPUT_DIR):
    raise Exception(f"Directory not found: {INPUT_DIR}")

files = sorted(glob.glob(os.path.join(INPUT_DIR, "*.txt")))

if len(files) == 0:
    raise Exception(
        f"No .txt files found in '{INPUT_DIR}'"
    )

# --------------------------------------------------
# Read data
# --------------------------------------------------

species_data = {}

for f in files:

    species = os.path.splitext(os.path.basename(f))[0]

    try:
        df = pd.read_csv(
            f,
            sep=r"\s+",
            header=None,
            names=["scaffold", "start", "end", "density"]
        )

        values = pd.to_numeric(
            df["density"],
            errors="coerce"
        ).dropna()

        if len(values) > 1:
            species_data[species] = values.values

    except Exception as e:
        print(f"Skipping {f}: {e}")

if len(species_data) == 0:
    raise Exception("No valid density data detected.")

# --------------------------------------------------
# Order species by median density
# --------------------------------------------------

ordered_species = sorted(
    species_data.keys(),
    key=lambda x: np.median(species_data[x]),
    reverse=True
)

# --------------------------------------------------
# Global x range
# --------------------------------------------------

all_values = np.concatenate(
    [species_data[s] for s in ordered_species]
)

xmin = np.min(all_values)
xmax = np.max(all_values)

padding = (xmax - xmin) * 0.05

xmin -= padding
xmax += padding

x = np.linspace(xmin, xmax, 1000)

# --------------------------------------------------
# Plot
# --------------------------------------------------

fig_height = max(5, len(ordered_species) * 0.45)

fig, ax = plt.subplots(
    figsize=(10, fig_height)
)

spacing = 1.0

for i, species in enumerate(ordered_species):

    vals = species_data[species]

    kde = gaussian_kde(vals)
    y = kde(x)

    y = y / y.max() * 0.8

    baseline = i * spacing

    # filled ridge
    ax.fill_between(
        x,
        baseline,
        baseline + y,
        alpha=0.8,
        linewidth=0
    )

    # outline
    ax.plot(
        x,
        baseline + y,
        linewidth=1
    )

    # median
    med = np.median(vals)

    ax.vlines(
        med,
        baseline,
        baseline + 0.75,
        linewidth=1.5
    )

# --------------------------------------------------
# Formatting
# --------------------------------------------------

ax.set_yticks(
    [i * spacing for i in range(len(ordered_species))]
)

ax.set_yticklabels(
    ordered_species,
    fontsize=10
)

ax.set_xlabel(
    "Gene density (genes per 100 kb)",
    fontsize=12
)

ax.set_ylabel("Species", fontsize=12)

ax.set_title(
    "Gene Density Distributions Across Genomes",
    fontsize=14
)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()

plt.savefig(
    "gene_density_ridgeline.png",
    dpi=600,
    bbox_inches="tight"
)

plt.savefig(
    "gene_density_ridgeline.pdf",
    bbox_inches="tight"
)

plt.show()

print(
    f"Ridgeline plot generated for {len(ordered_species)} species."
)
