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

   


   
