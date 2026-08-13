# TE Divergence visualization

This script can be used for visualizing `transposable elements (TE) divergence` statistics across species. The script takes a tab-delimited summary file containing TE divergence statistics and generates visualizations. The figures are exported in both PDF and SVG formats.
The script is designed as a visualization tool following the [`1.repeat-landscapes.py`](./1.repeat-landscapes.py), using the output `TE_divergence_summary.tsv`.

## Requirements
Python 3.8 or newer.

Required Python packages;

`pandas` `numpy` `matplotlib`

To install;
```
pip install pandas numpy matplotlib
```

## Usage
```
python3 2.plot_TE-divergence-boxplots.py -i TE_divergence_summary.tsv
```
`-i` / `--input` TE divergence summary TSV file

## Input file

### 1. Tab-delimited transposable elements divergence

This file is one of the output generated from previous repeat landscapes mapping script [`1.repeat-landscapes.py`](./1.repeat-landscapes.py)
