# Hox gene identification and phylogenetic analysis

This script is used to identify putative Hox genes from multiple peptide FASTA files, combining candidate sequences with reference Hox proteins, performing multiple sequence alignment and removing sequences that does not contain informative amino-acid residues after alignment. The program also reconstructs a phylogenetic tree using IQ-TREE.


## Requirements
Python 3.8 or newer
External programs required; `BLAST+` `MAFFT` `IQTREE2/3`

## Conda environment

Create a new conda environment:

`conda create -n hoxpull python=3.10` conda activate hoxpull

Install `BLAST+` `MAFFT` and `IQ-TREE`

```
conda install -c bioconda -c conda-forge blast mafft iqtree
```
Verify the versions of BLAST, MAFFT and IQ-TREE and also the path to IQ-TREE as the script requires the path to IQ-TREE.

`which iqtree` will show the path for the same

## Usage
The script accepts six command line arguments.

```
python3 1.Hox-multi-extract.py -p .pep_directory -r hox_references -i path_to_iq-tree -o outdir -t threads -e e-value -m max-targets
```
`-p` path to directory containing all peptides
`-r` reference hox peptides .fasta
`-i` path to IQ-TREE executable
`-o` output directory
`-t` number of CPU threads
`-e` BLASTP E-value threshold
`-m` Maximum blast target sequences


# Input files

## 1. 


 
