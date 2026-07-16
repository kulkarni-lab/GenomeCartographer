#!/usr/bin/env python3

"""
Ridgeline (joyplot) of gene density distributions across genomes.
Outputs metrics (raw probability density) with localized scale bars.
"""

import os
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

# --------------------------------------------------
# Input configuration
# --------------------------------------------------
INPUT_DIR = "gene_density"
if not os.path.isdir(INPUT_DIR):
    raise Exception(f"Directory not found: {INPUT_DIR}")

files = sorted(glob.glob(os.path.join(INPUT_DIR, "*.txt")))
if not files:
    raise Exception(f"No .txt files found in '{INPUT_DIR}'")

# --------------------------------------------------
# Read and process data
# --------------------------------------------------
species_data = {}
for f in files:
    species = os.path.splitext(os.path.basename(f))[0]
    try:
        df = pd.read_csv(f, sep=r"\s+", header=None, names=["scaffold", "start", "end", "density"])
        values = pd.to_numeric(df["density"], errors="coerce").dropna()
        if len(values) > 1:
            species_data[species] = values.values
    except Exception as e:
        print(f"Skipping {f}: {e}")

if not species_data:
    raise Exception("No valid density data detected.")

# Order by median density
ordered_species = sorted(species_data.keys(), key=lambda x: np.median(species_data[x]), reverse=True)

# --------------------------------------------------
# Setup Global Ranges & Find Max Density
# --------------------------------------------------
all_values = np.concatenate([species_data[s] for s in ordered_species])
x = np.linspace(all_values.min(), all_values.max(), 1000)

# Calculate global peak height to determine spacing and scale bars
max_density = 0
for species in ordered_species:
    kde = gaussian_kde(species_data[species])
    max_density = max(max_density, np.max(kde(x)))

# Output the maximum density to the terminal as requested
print(f"Maximum calculated probability density: {max_density}")

# --------------------------------------------------
# Plotting Setup
# --------------------------------------------------
fig_height = max(6, len(ordered_species) * 0.7)
fig, ax = plt.subplots(figsize=(10, fig_height))

# Spacing between baselines
spacing = max_density * 0.8 

# Scale bar configuration (e.g., 50% of the max density height, formatted nicely)
# Using 2 significant digits for a clean label
scale_height = float(f"{max_density * 0.5:.2g}") 

# X-axis configuration for scale bars (placed slightly left of the data minimum)
x_range = x.max() - x.min()
scale_x = x.min() - (x_range * 0.03)
tick_w = x_range * 0.01  # Width of the horizontal caps on the scale bar

# Extend the x-axis limit so the scale bars have room and don't overlap y-axis text
ax.set_xlim(x.min() - (x_range * 0.1), x.max() + (x_range * 0.05))

# --------------------------------------------------
# Draw Ridgelines
# --------------------------------------------------
for i, species in enumerate(ordered_species):
    vals = species_data[species]
    kde = gaussian_kde(vals)
    y = kde(x)
    baseline = i * spacing

    # Filled ridge
    ax.fill_between(x, baseline, baseline + y, alpha=0.7, linewidth=0.5, edgecolor="black")

    # Median indicator
    med = np.median(vals)
    med_y = kde(med)[0]
    ax.vlines(med, baseline, baseline + med_y, color='red', linestyle='--', linewidth=1.5)

    # --- Localized Scale Bar ---
    # Vertical line
    ax.plot([scale_x, scale_x], [baseline, baseline + scale_height], color='black', linewidth=1)
    # Bottom cap
    ax.plot([scale_x - tick_w, scale_x], [baseline, baseline], color='black', linewidth=1)
    # Top cap
    ax.plot([scale_x - tick_w, scale_x], [baseline + scale_height, baseline + scale_height], color='black', linewidth=1)
    
    # Scale bar text label
    ax.text(
        scale_x - tick_w * 1.5, 
        baseline + (scale_height / 2), 
        f"{scale_height}", 
        va='center', 
        ha='right', 
        fontsize=9,
        color='black'
    )

# --------------------------------------------------
# Formatting
# --------------------------------------------------
# Y-ticks represent the baselines for the species names
ax.set_yticks([i * spacing for i in range(len(ordered_species))])
ax.set_yticklabels(ordered_species, fontsize=10)

ax.set_xlabel("Gene density (genes per 100 kb)", fontsize=12)
# Removed standard y-axis label since the scale bars now explicitly state the metric
ax.set_title("Gene Density Distributions (with Probability Density Scales)", fontsize=14)

# Hide unnecessary spines for a clean look
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_visible(False) # Left spine hidden to emphasize individual scale bars

# Ensure ticks on the bottom are preserved but hidden on the left
ax.tick_params(axis='y', length=0) 

plt.tight_layout()
plt.savefig("gene_density_ridgeline_metric.png", dpi=600, bbox_inches="tight")
plt.savefig("gene_density_ridgeline_metric.pdf", bbox_inches="tight")

plt.show()
