# Genome architecture analysis

This script is a pipeline for comparative analysis of genome architecture using `GTF`, `GFF` and `GFF3` gene annotations.

The script extracts gene, exon and intron structures from genome annotations and compares gene architecture across multiple speciesa and generates `violin` and `boxplot` visualization.

The script calculates;
```
Exon length
Intron length
Gene span
Number of exons per gene
Number of introns per gene
Total exon length
Gene compactness
Species-level summary statistics
```

## Requirements
Python 3.8 or newer is recommended.

Required Python packages;
`numpy` `pandas` `matplotlib` `scipy` `statsmodels`

To install;
```
pip install numpy pandas matplotlib scipy statsmodels
```


## Usage
```
python3 1.Intron-from-gtf_corr.py annotations genome_architecture order.txt
```
`annotations`         path to directory containing GTF/GFF/GFF3 files

`genome architecture` prefix used for all output files

`order.txt`           path to text file containing species names in the desired order

## Input files

## 1. Genome annotation files
   The script accepts standard 9-column; `.gtf` `.gff` `.gff3`

   The script uses following feature types; `gene` `mRNA` `transcript` `exon` `CDS`
   
   All annotations files should be placed in a single directory.

The directory ./annotations should contain files named with one of the above extensions.

```
Example: ./annotations/
          species1.gff3
          species2.gff3          
          species3.gff3
```

The filename determines the species name.

The parser attempts to identify transcript and gene IDs using several common attributes:
```
transcript_id
gene_id
Parent
ID
gene
```

Example annotaion:

```
scaffold1	source	gene	1000	10000	.	+	.	gene_id "gene001"
scaffold1	source	mRNA	1000	10000	.	+	.	gene_id "gene001"; transcript_id "tx001"
scaffold1	source	exon	1000	1500	.	+	.	gene_id "gene001"; transcript_id "tx001"
scaffold1	source	exon	3000	3500	.	+	.	gene_id "gene001"; transcript_id "tx001"
scaffold1	source	exon	7000	10000	.	+	.	gene_id "gene001"; transcript_id "tx001"
```

## 2. Species order file
   The TXT file defines the species order and reference species

Example: [`order.txt`](Intron_Lengths/order.txt)
```
Indi2
Ixod
Ixpers
Rhip
Avul
Otur
Nlon
Phytos
Galen
Varr
```

The first species is treated as the reference and the order also controls the ordering of species in the output plots. The order file should contain one species name per line.


# Output
The script generates five tab-delimited tables and three SVG files followed by the prefix given in the command as filename;

```
genome_architecture_introns.tsv
genome_architecture_exons.tsv
genome architecture_compactness.tsv
genome_architecture_gene_summary.tsv
genome_architecture_species_summary.tsv
genome_architecture_intron_violin.svg
genome_architecture_exon_boxplot.svg
genome_architecture_compactness_boxplot.svg
```

### Output files description

#### TABLES

`genome_architecture_introns.tsv`

- Individual intron measurements: `species` `gene` `intron_length`

`genome_architecture_exons.tsv`

- Individual exon measurements: `species` `gene` `exon_length`

`genome architecture_compactness.tsv`

- Gene-level compactness: `species` `gene` `compactness`

`genome_architecture_gene_summary.tsv`

- Complete gene-level summary: `species` `gene` `gene_span` `num_exons` `num_introns` `total_exon_length` `compactness`

`genome_architecture_species_summary.tsv`

- Species-level summary: `species` `median_intron_length` `mean_intron_length` `median_exon_length` `mean_exon_length` `median_compactness` `mean_compactness` `mean_introns_per_gene` `median_gene_span` `num_genes`

#### FIGURES

`genome_architecture_intron_violin.svg` - **Intron size distribution**

- A violin plot showing the distribution of `intron lengths` across species.
- The y-axis uses a logarithmic scale.

`genome_architecture_exon_boxplot.svg` - **Exon length distribution**

- A boxplot showing `exon-length` distributions.
- The y-axis uses a logarithmic scale.

`genome_architecture_compactness_boxplot.svg` - **Gene compactness**

- A boxplot showing `gene compactness` across species.
- The y-axis is constrained to approximately: `0–1`

## Calculations and input selections

### 1. Transcript Selection
Many genes have multiple transcript isoforms. To avoid treating every transcript as an independent gene, the script selects one transcript per gene. For each gene, script calculates the total `exon/CDS` length for every transcript and retains the transcript with the largest total feature length.

Example:
```
Gene_01

transcript1 = 1,200 bp
transcript2 = 2,500 bp
transcript3 = 1,900 bp

The script retains:
transcript2 for analysis
```

This produces one representative transcript per gene.

### 2. Gene Span

Gene span is calculated as:
`gene span = maximum feature end - minimum feature start + 1`

```
For example:

minimum coordinate = 100,000
maximum coordinate = 150,000

gives:

gene span = 50,001 bp
```

### 3. Exon Length

Each exon is measured as:
`exon length = end - start + 1`

The resulting measurements are stored in: `*_exons.tsv`

### 4.Intron Length
Introns are inferred from the gaps between consecutive exon features.

For two adjacent exons: 
```
exon 1: 1000–1500
exon 2: 3000–3500

the inferred intron length is: 1501–2999 = 1499 bp

```
Only positive-length gaps are recorded as introns. Inferred introns are stored in:
`*_introns.tsv`


### 5. Compactness
Gene compactness is calculated as:

`compactness = total exon length / gene span`
```
For example:

total exon length = 10,000 bp
gene span = 50,000 bp
compactness = 0.20
```
Higher values indicate that a larger fraction of the gene span consists of exonic sequence. Lower values indicate a greater contribution of intronic sequence to the gene span.

## Statistical Analysis
The script compares every species against the reference species defined by the first line of `order.txt`. For each metric, the script uses the Wilcoxon rank-sum test: `scipy.stats.ranksums`

The comparisons are:
```
species2 vs reference
species3 vs reference
species4 vs reference
...
The resulting P-values are corrected using Bonferroni correction:

multipletests(
    ps,
    method="bonferroni"
)
```
Only adjusted P-values below: `0.05` are displayed on the figures.

## Significance Labels

The figures use:
```
*       p < 0.05
**      p < 0.01
***     p < 0.001
****    p < 0.0001
```
These labels are based on the Bonferroni-adjusted P-values.
