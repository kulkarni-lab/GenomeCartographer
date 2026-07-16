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

    print(f"Reading {f}")

    # Skip metadata lines beginning with #
    df = pd.read_csv(
        f,
        comment="#"
    )

    species = df["species"].iloc[0]

    # Keep only family + filtered hits
    df = df[["family", "filtered_hits"]]

    # Convert to numeric copy numbers
    df["filtered_hits"] = (
        pd.to_numeric(
            df["filtered_hits"],
            errors="coerce"
        )
        .fillna(0)
        .astype(int)
    )

    # Rename column to species name
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
# CUSTOM COLOR MAP
###########################################################
# 0 = white
# 1-5 = increasing red intensity
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

fig, ax = plt.subplots(
    figsize=(10,14)
)

im = ax.imshow(
    plot_df.values,
    aspect='auto',
    interpolation='none',   # IMPORTANT FOR VECTOR CELLS
    cmap=cmap,
    norm=norm
)

###########################################################
# DRAW CELL BORDERS
###########################################################

ax.set_xticks(
    [x - 0.5 for x in range(1, len(plot_df.columns))],
    minor=True
)

ax.set_yticks(
    [y - 0.5 for y in range(1, len(plot_df.index))],
    minor=True
)

ax.grid(
    which="minor",
    color="grey",
    linestyle='-',
    linewidth=0.3
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

cbar.set_label(
    "miRNA copy number",
    fontsize=10
)

###########################################################
# TITLE
###########################################################

ax.set_title(
    "miRNA Copy Number Evolution Across Acari",
    fontsize=14
)

plt.tight_layout()

###########################################################
# SAVE COPY NUMBER FIGURES
###########################################################

plt.savefig(
    "Acari_miRNA_copy_number_heatmap.pdf",
    bbox_inches='tight'
)

plt.savefig(
    "Acari_miRNA_copy_number_heatmap.svg",
    bbox_inches='tight'
)

plt.savefig(
    "Acari_miRNA_copy_number_heatmap.png",
    dpi=400,
    bbox_inches='tight'
)

plt.close()

###########################################################
# PRESENCE / ABSENCE MATRIX
###########################################################

binary_df = (heatmap_df > 0).astype(int)

###########################################################
# SORT miRNAs BY CONSERVATION
###########################################################

binary_df["conservation"] = binary_df.sum(axis=1)

binary_df = binary_df.sort_values(
    by="conservation",
    ascending=False
)

plot_df2 = binary_df.drop(
    columns=["conservation"]
)

###########################################################
# CONSERVATION HEATMAP
###########################################################

fig, ax = plt.subplots(
    figsize=(10,14)
)

im = ax.imshow(
    plot_df2.values,
    aspect='auto',
    interpolation='none',
    cmap=ListedColormap([
        "#ffffff",
        "#08519c"
    ]),
    vmin=0,
    vmax=1
)

###########################################################
# CELL BORDERS
###########################################################

ax.set_xticks(
    [x - 0.5 for x in range(1, len(plot_df2.columns))],
    minor=True
)

ax.set_yticks(
    [y - 0.5 for y in range(1, len(plot_df2.index))],
    minor=True
)

ax.grid(
    which="minor",
    color="grey",
    linestyle='-',
    linewidth=0.3
)

###########################################################
# LABELS
###########################################################

ax.set_xticks(range(len(plot_df2.columns)))

ax.set_xticklabels(
    plot_df2.columns,
    rotation=45,
    ha='right',
    fontsize=9
)

ax.set_yticks(range(len(plot_df2.index)))

ax.set_yticklabels(
    plot_df2.index,
    fontsize=7
)

###########################################################
# TITLE
###########################################################

ax.set_title(
    "miRNA Family Conservation Across Acari",
    fontsize=14
)

plt.tight_layout()

###########################################################
# SAVE CONSERVATION FIGURES
###########################################################

plt.savefig(
    "Acari_miRNA_conservation_heatmap.pdf",
    bbox_inches='tight'
)

plt.savefig(
    "Acari_miRNA_conservation_heatmap.svg",
    bbox_inches='tight'
)

plt.savefig(
    "Acari_miRNA_conservation_heatmap.png",
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

print("\nCopy number heatmap:")
print(" - Acari_miRNA_copy_number_heatmap.pdf")
print(" - Acari_miRNA_copy_number_heatmap.svg")
print(" - Acari_miRNA_copy_number_heatmap.png")

print("\nConservation heatmap:")
print(" - Acari_miRNA_conservation_heatmap.pdf")
print(" - Acari_miRNA_conservation_heatmap.svg")
print(" - Acari_miRNA_conservation_heatmap.png")
