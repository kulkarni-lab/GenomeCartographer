# Hox gene cluster visualization

This script is used to extract genomic coordinates of Hox genes from chromosome/scaffold coordinate files and generate high-quality diagrams for **Hox genes** across multiple genomes.

This script is designed as a visualization tool following the [`1.Hox-multi-extract.py`](../Hox-extraction/1.Hox-multi-extract.py) script.


## Requirements
Python 3.8 or newer is recommended.

Required Python packages;
`pandas` `matplotlib`

To install;
```
pip install pandas matplotlib
```


## Usage
```
python3 3.hox-plot.py -c input.csv -genomedir ./genomes -chromdir ./chroms -o output-directory
```
`-c` a comma separated .csv input file with four column
`-genome` path to directory where all genome files are stored
`-chromdir` path to directory where all chrom files are stored
`-o` output directory name

# Input files

**The script requires three types of input files**

## 1. Hox gene CSV file
Example: `hox_genes.csv`

- This file is generated using the output files created from [1.Hox-multi-extract.py](../Hox-extraction/1.Hox-multi-extract.py).
- The treefile (hox.treefile) needs to be manually examined to find monophyletic gene groups that cluster together with the reference fasta sequence.
- It must contain four columns: `genome_name,chrom_file,sequence_name,hoxgene`

`genome_name` - Genome FASTA filename
`chrom_file` - Corresponding chromosome/scaffold coordinate file
`sequence_name` - Gene ID for the Hox gene identified
`hoxgene` - Hox gene name

Example:
```
genome_name,chrom_file,sequence_name,hoxgene
Species_A.fna,Species_A.chrom,g1.t1,lab
Species_A.fna,Species_A.chrom,g55.t1,pb
Species_A.fna,Species_A.chrom,g212.t1,abdA
Species_A.fna,Species_A.chrom,g310.t1,scr
Species_A.fna,Species_A.chrom,g450.t1,antp
```

## 2. Genome FASTA directory
- Argument - `-genomedir`

- The directory should contain all the genome FASTA files listed in `hox_genes.csv`.
- The genome filename should match exactly with the name in `hox_genes.csv`.


## 3. Chromosome coordinate files
- Argument - `-chromdir`

- The directory should contain all chromosome/scaffold coordinate files.
- each `.chrom` file should have five whitespace-seperated columns: `gene_ID` `scaffold` `strand` `start` `end`

```
Example:
XP_054165658.1  NW_026533251.1  -       43823   45701
XP_054159110.1  NW_026533221.1  +       1267302 1269922
XP_054167116.1  NW_026533194.1  +       314816  317216
XP_054154760.1  NW_026533201.1  -       1159699 1161490
XP_054168353.1  NW_026533196.1  +       686881  689953
XP_054158519.1  NW_026533191.1  -       1257011 1257776
XP_054163170.1  NW_026533243.1  +       2120421 2121405
XP_054167320.1  NW_026533253.1  -       1477556 1478522
XP_054156547.1  NW_026533210.1  +       3481755 3483330
XP_054160724.1  NW_026533221.1  -       6296025 6302065
```

- The script automatically simplifies scaffold names.
- Example: `CM075101.1_Phytoseiulus_persimilis_chromosome_1` becomes `CM075101.1`


# Output

The script generates three main files;
```
hox_clusters.annotated.csv
hox_cluster.png
hox_cluster.pdf
```

### Output files description

**1** `hox_clusters.annotated.csv`

- This CSV file contains all the coordinates actually used for plotting.
- The CSV file contains:
```
genome
chrom_file
sequence
hoxgene
scaffold
start
end
strand
```

**2** `hox_clusters.png`
**3** `hox_clusters.pdf`

**2** and **3** are graphical files both PNG and PDF format respectively. PNG file is 300 DPI and PDF file is vector-based, which can be used for easier editing.

Example:

<img width="1010" height="560" alt="hox_clusters" src="https://github.com/user-attachments/assets/d616d1a1-f99c-4f5a-a305-05ed414488b9" />





