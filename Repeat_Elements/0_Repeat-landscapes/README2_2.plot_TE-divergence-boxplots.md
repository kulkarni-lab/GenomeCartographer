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

- This file is one of the output generated from previous repeat landscapes mapping script [`1.repeat-landscapes.py`](./1.repeat-landscapes.py)
- It should contain following columns;

```
Species
Mean_divergence
Median_divergence
Q1_divergence
Q3_divergence
IQR_divergence
SD_divergence
```

Example: 
```
Species	Mean_divergence	Median_divergence	Q1_divergence	Q3_divergence	IQR_divergence	SD_divergence
speciesA	12.43	10.00	4.00	18.00	14.00	9.21
speciesB	15.21	13.00	7.00	21.00	14.00	10.32
```

# Output

- The script generates eight files: four PDFs and four PNG files.

```
Median_IQR_divergence.pdf
Median_IQR_divergence.svg
TE_divergence_box_summary.pdf
TE_divergence_box_summary.svg
TE_divergence_bubbleplot.pdf
TE_divergence_bubbleplot.svg
TE_divergence_heatmap.pdf
TE_divergence_heatmap.svg
```
### Output files description

#### 1. Median divergence ± IQR
```
Median_IQR_divergence.pdf
Median_IQR_divergence.svg
```
This figure displays: 
- Median TE divergence as the central point.
- Q1 and Q3 as the lower and upper limits of the error bar

#### 2. Box-style summary
```
TE_divergence_box_summary.pdf
TE_divergence_box_summary.svg
```
This visualization represents:
- Rectangle → Q1 to Q3
- Horizontal line → Median
- Vertical line → Median ± SD

#### 3. TE age structure bubble plot
```
TE_divergence_bubbleplot.pdf
TE_divergence_bubbleplot.svg
```
This plot compares:
- X-axis → Median TE divergence
- Y-axis → Mean TE divergence
- Bubble size → IQR

`X = Median divergence ` `Y = Mean divergence` `Size = IQR`

**Interpretation**

This visualization allows comparison of both the central tendency and variability of TE divergence.

For example:
```
High median + high mean → relatively high divergence
Low median + low mean → relatively low divergence
Large bubble → greater variability
Small bubble → narrower divergence distribution
```

#### 4. Summary statistics heatmap
```
TE_divergence_heatmap.pdf
TE_divergence_heatmap.svg
```

The heatmap displays: 
```
 Mean divergence
 Median divergence
 Q1 divergence
 Q3 divergence
 IQR divergence
 SD divergence
```
- Heatmap provides a quick overview of how the different summary statistics vary across species.
