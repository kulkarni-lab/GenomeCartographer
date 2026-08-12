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
