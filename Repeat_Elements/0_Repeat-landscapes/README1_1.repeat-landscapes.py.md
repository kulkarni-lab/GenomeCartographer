# RepeatMasker TE landscape analysis

This script helps in generating transposable elements (TE) landscapes from RepeatMasker `outputfiles` across multiple genomes. The pipeline calculates `repeat divergence statistics`, estimates the percentage of each genome occupied by `repetitive DNA` across divergence bins, generates individual `TE-class landscapes` for each species, and produces a combined `multi-species repeat landscape`.

## Requirements 
Python 3.8 or newer.

Required python packages;
`pandas` `numpy` `matplotlib` `scipy`

To install;
```
pip install pandas numpy matplotlib scipy
```


## Usage
```
python3 1.repeat-landscapes.py -i input -g genome_size_file
```
`-i` / `--indir` path to RepeatMasker `.out` files

`-g` / `--genomes` species genome size file table

# Input files

## 1. RepeatMasker output directory
   Directory should be supplied using ;
   `-i` or `--indir`

Place all RepeatMasker `.out` files in one directory.
```
Example: ./repeatmasker.out/
         species1.out
         species2.out
         species3.out
         species4.out
```

The filename determines the species name.

## 2. Genome size table
   The genome-size file must contain at least two columns: `species` `genome_size`. It must be tab separated file.

   Example: [`gen-sizes.txt`](./gen-sizes.txt)
```
species   genome_size
MESOSTIG_Neoseiulus.fna   199210599
MESOSTIG_Phytoseiulus_persimilis.fna   181015949
Galend.fna   151460383
ACARI_Argas_vulgaris.fna   1340163888
```

- Genome size must be given in **base pairs**
- The species name must be matching exactly with the `.out` filename and `genome-size table`.  
[!NOTE] files without matching names are skipped.

# Output
The script generates two figures (PNG and PDF) for each successfully processed species and all species TE landscape as well.

```
speciesA_landscape_percentGenome.pdf
speciesA_landscape_percentGenome.svg
speciesB_landscape_percentGenome.pdf
speciesB_landscape_percentGenome.svg
speciesC_landscape_percentGenome.pdf
speciesC_landscape_percentGenome.svg
All_species_TE_landscape_percentGenome.pdf
All_species_TE_landscape_percentGenome.svg
TE_divergence_summary.tsv
```

### Output files description

`speciesA_landscape_percentGenome.pdf`
`speciesA_landscape_percentGenome.svg`

- Individual species plots detailed TE curves for `DNA` `LINE` `LTR` `SINE` `Helitron` `Unknown`.
- The `x-axis` is percent divergence from consensus and `y-axis` is percent genome occupied.


`All_species_TE_landscape_percentGenome.pdf`
`All_species_TE_landscape_percentGenome.svg`

- This figure shows the overall repeat landscapes for all species.
- Combined plot does not separate curves by TE class, instead represents the total repeat occupancy across divergence bins.

`TE_divergence_summary.tsv`

The script generates a summary table with the following columns;
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

## How the script works

#### RepeatMasker _.out_ Parsing

- The script extracts the following information from each **RepeatMasker** record: `divergence from consensus` `genomic start` `genomic end``repeat classification`

**Repeat length is calculated as:** `end - start + 1`

[!NOTE] All the **RepeatMasker** header lines and malformed records are ignored.

### Divergence Analysis

- **RepeatMasker** reports divergence between individual repeat copies and their consensus sequence. Repeat length is used as the weighting factor.

- The script calculates several statistics:
```
weighted mean divergence
weighted median divergence
weighted first quartile
weighted third quartile
weighted interquartile range
weighted standard deviation
```

### Weighted mean

- The weighted mean is approximately:
```
Σ(divergence × repeat length)
--------------------------------
       Σ(repeat length)
```
- This means long repeat sequences contribute more strongly than short repeat sequences.

### Weighted Median

- The `median divergence` is also calculated using repeat length as the weight. This represents the divergence value at which approximately 50% of the total repeat sequence lies below and above that value.
- The same approach is used to calculate:
```
Q1 = 25th percentile
Q3 = 75th percentile
and:
IQR = Q3 - Q1
```
### Divergence Binning

- The repeat landscape is calculated using 1% divergence bins: `0–1%` `1–2%` `2–3%` `...` `39–40%`

- The x-coordinate is represented by the bin midpoint: `0.5` `1.5` `2.5` `...` `39.5`

Therefore, the x-axis represents: Percent divergence from consensus

### Genome Occupancy

- For each divergence bin, the script calculates:
```
repeat base pairs in bin
------------------------ × 100
genome size
````
The resulting value is: `% genome occupied`

Example: 
```
If;
repeat sequence in bin = 1,000,000 bp
genome size = 100,000,000 bp

then: genome occupancy = 1%
```
### TE Landscape

- Each species receives an individual `TE landscape`. The individual landscape contains curves for:
```
DNA
LINE
LTR
SINE
Helitron
Unknown
```
The x-axis is: `Percent divergence from consensus` ,the y-axis is: `% genome occupied`

### Smoothing

- The landscape curves are smoothed using a `Savitzky–Golay` filter: `savgol_filter()`. The smoothing is used for visualization.


with:

window = 7 bins

polynomial order = 3

- The underlying divergence-bin occupancy values remain the basis for the peak calculation.

### Peak Divergence

- For each species, the script identifies the divergence bin with the highest total repeat occupancy.

Example:

Peak divergence   : 4.5%

Peak occupancy    : 2.3814% genome

Which means the largest amount of repeat sequence occurs in the: `4–5%` divergence interval. Peak detection is performed on the unsmoothed binned values.

## Interpreting TE Landscapes

- The general interpretation is:

`Low divergence`

- Repeats at low divergence from consensus are generally more similar to their inferred consensus sequence. 
- A strong peak at low divergence may indicate a substantial contribution from relatively young or recently amplified repeat copies.

`High divergence`

- Repeats at higher divergence have accumulated more substitutions relative to their consensus.
- A broad distribution extending toward high divergence indicates a substantial amount of older or highly diverged repetitive sequence.
