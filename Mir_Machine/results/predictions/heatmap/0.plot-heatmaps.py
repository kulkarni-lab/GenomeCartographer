import pandas as pd
import glob
from functools import reduce
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm

###########################################################
# OPTIONAL TAXON ORDER FILE
###########################################################
# taxon_order.txt should contain ONE species per line
#
# Example:
#
# Opilioacarus
# Ixodes
# Derm
# Mesostigma
#
###########################################################

USE_TAXON_ORDER = True
TAXON_ORDER_FILE = "taxon_order.txt"

###########################################################
# READ TAXON ORDER
###########################################################

taxon_order = None

if USE_TAXON_ORDER:

    with open(TAXON_ORDER_FILE) as infile:

        taxon_order = [
            line.strip()
            for line in infile
            if line.strip()
        ]

###########################################################
# READ ALL MIRMachine FILES
###########################################################

files = glob.glob("*.heatmap.csv")

dfs = []

for f in files:

    df = pd.read_csv(
        f,
        comment="#"
    )

    species = df["species"].iloc[0]

    df = df[["family", "filtered_hits"]]

    df["filtered_hits"] = (
        pd.to_numeric(
            df["filtered_hits"],
            errors="coerce"
        )
        .fillna(0)
        .astype(int)
    )

    df = df.rename(
        columns={"filtered_hits": species}
    )

    dfs.append(df)

###########################################################
# MERGE ALL SPECIES
###########################################################

merged = reduce(
    lambda left, right:
    pd.merge(
        left,
        right,
        on="family",
        how="outer"
    ),
    dfs
)

merged = merged.fillna(0)

###########################################################
# ENSURE NUMERIC
###########################################################

for col in merged.columns[1:]:

    merged[col] = (
        pd.to_numeric(
            merged[col],
            errors="coerce"
        )
        .fillna(0)
        .astype(int)
    )

###########################################################
# SAVE RAW MATRIX
###########################################################

merged.to_csv(
    "combined_miRNA_copy_numbers.csv",
    index=False
)

###########################################################
# PREPARE MATRIX
###########################################################

heatmap_df = merged.set_index("family")

###########################################################
# REORDER TAXA IF PROVIDED
###########################################################

if taxon_order:

    # Keep only taxa present in matrix
    valid_taxa = [
        t for t in taxon_order
        if t in heatmap_df.columns
    ]

    heatmap_df = heatmap_df[valid_taxa]

###########################################################
# CAP VALUES AT 5+
###########################################################

plot_df = heatmap_df.copy()

plot_df[plot_df >= 5] = 5

###########################################################
# COLOR MAP
###########################################################

colors = [
    "#ffffff",  # 0
    "#fee5d9",  # 1
    "#fcae91",  # 2
    "#fb6a4a",  # 3
    "#de2d26",  # 4
    "#a50f15"   # 5+
]

cmap = ListedColormap(colors)

bounds = [0,1,2,3,4,5,6]

norm = BoundaryNorm(bounds, cmap.N)

###########################################################
# COPY NUMBER HEATMAP
###########################################################

fig, ax = plt.subplots(figsize=(10,14))

im = ax.imshow(
    plot_df.values,
    aspect='auto',
    interpolation='nearest',
    cmap=cmap,
    norm=norm
)

###########################################################
# LABELS
###########################################################

ax.set_xticks(range(len(plot_df.columns)))

ax.set_xticklabels(
    plot_df.columns,
    rotation=45,
    ha='right',
    fontsize=9
)

ax.set_yticks(range(len(plot_df.index)))

ax.set_yticklabels(
    plot_df.index,
    fontsize=7
)

###########################################################
# COLORBAR
###########################################################

cbar = plt.colorbar(
    im,
    ax=ax,
    ticks=[0.5,1.5,2.5,3.5,4.5,5.5]
)

cbar.ax.set_yticklabels([
    "0",
    "1",
    "2",
    "3",
    "4",
    "5+"
])

cbar.set_label("miRNA copy number")

###########################################################
# TITLE
###########################################################

ax.set_title(
    "miRNA Copy Number Evolution Across Acari",
    fontsize=14
)

plt.tight_layout()

###########################################################
# SAVE COPY NUMBER FIGURE
###########################################################

plt.savefig(
    "Acari_miRNA_copy_number_heatmap.pdf",
    bbox_inches='tight'
)

plt.savefig(
    "Acari_miRNA_copy_number_heatmap.png",
    dpi=400,
    bbox_inches='tight'
)

plt.close()

###########################################################
# CONSERVED miRNAs
###########################################################

binary_df = (heatmap_df > 0).astype(int)

conserved = binary_df[
    binary_df.sum(axis=1) == len(binary_df.columns)
]

###########################################################
# CONSERVED FIGURE
###########################################################

fig, ax = plt.subplots(
    figsize=(8, max(4, len(conserved)*0.25))
)

im = ax.imshow(
    conserved.values,
    aspect='auto',
    interpolation='nearest',
    cmap=ListedColormap(["white", "#08519c"]),
    vmin=0,
    vmax=1
)

###########################################################
# LABELS
###########################################################

ax.set_xticks(range(len(conserved.columns)))

ax.set_xticklabels(
    conserved.columns,
    rotation=45,
    ha='right',
    fontsize=9
)

ax.set_yticks(range(len(conserved.index)))

ax.set_yticklabels(
    conserved.index,
    fontsize=8
)

###########################################################
# TITLE
###########################################################

ax.set_title(
    "Conserved miRNA Families Across Acari",
    fontsize=14
)

plt.tight_layout()

###########################################################
# SAVE CONSERVED FIGURE
###########################################################

plt.savefig(
    "Acari_conserved_miRNAs.pdf",
    bbox_inches='tight'
)

plt.savefig(
    "Acari_conserved_miRNAs.png",
    dpi=400,
    bbox_inches='tight'
)

plt.close()

###########################################################
# SUMMARY
###########################################################

print("\nFinished!\n")

print("Generated files:")

print(" - combined_miRNA_copy_numbers.csv")
print(" - Acari_miRNA_copy_number_heatmap.pdf")
print(" - Acari_miRNA_copy_number_heatmap.png")
print(" - Acari_conserved_miRNAs.pdf")
print(" - Acari_conserved_miRNAs.png")
