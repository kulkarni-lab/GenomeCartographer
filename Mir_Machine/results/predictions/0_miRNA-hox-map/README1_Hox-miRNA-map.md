# Hox-miRNA mapping

This script maps miRNA loci to HOX gene clusters across multiple taxa. 
The script reads HOX gene coordinates and miRNA GFF annotations, determines
the genomic relationship between each miRNA and the corresponding HOX cluster, 
calculates an HOX association score, generates presence/absence matrices, extracts 
canonical HOX-associated miRNAs, and produces a heatmap of HOX association scores.

## Requirements
Python 3.8 or newer.

Required Python packages;
`pandas`
`numpy`
`matplotlib`

To install;
```
pip install pandas numpy matplotlib
```


## Usage
```
python3 1.hox-miRNA-map.py -m input_of_hox-miRNAs.gffs --hox hox_clusters.annotated.csv
```
`-m` path to gffs_hox_miRNA 

`--hox` path to hox coordinate CSV file

# Input files

## 1. miRNA GFF files 
   Directory should be supplied using;
   `-m` or `--mirna_dir`

The directory (gffs_hox-miRNA) should contain files named as:

`filename_miRNA.gff`

```
Example: Avul_miRNA.gff
         Galen_miRNA.gff
         Ixod_miRNA.gff
```        
  The script needs regular 9 column `gff/gff3`

  Example: [`gffs_hox-miRNA`](./gffs_hox-miRNA)
```
CM147321.1	MirMachine	microRNA	138807572	138807631	50.3	+	.	gene_id=Bantam.PRE;E-value=5.2e-05;sequence_with_30nt=TAACACGACTGGTGGAGGATCAGACAAAACTGGTTTTCACAATGATCATCCAGATGTGTCCGATATCTGAGATCATTGTGAAAGCTGATTTTGTTGTTTCGACAACGAGGGAAGCGGGAC;seed=(p3_seed(GAGATCA*))`
CM147319.1	MirMachine	microRNA	208675234	208675291	84.6	+	.	gene_id=Iab-4.PRE;E-value=7.3e-17;sequence_with_30nt=GACGGCGATGGCTGGATCGGCATCTCCTGTTCGTATACTGAGTGTATCCTGAGTGGACAACTTTCCGGTATACCTTCAGTATACGTAACAGGCGACCCGTTTCGAGAGCGACTTGGAT;seed=(p5_seed(CGTATAC*))
CM147320.1	MirMachine	microRNA	73825156	73825227	67.4	+	.	gene_id=Let-7.PRE;E-value=3e-10;sequence_with_30nt=TTAAATCATCGCTTTTGCGTGCGCTCAGTGTGAGGTAGTAGGTTGTATAGTTGAGAACTACATCTGTCTCGGAGGCCTAACTGTACAACTTGCTAACTTACTCTGCGTGGCACTGTTAAACGTCGCCTCTCG;seed=(p5_seed(TGAGGTA,GAGGTAG*))
```
## 2. HOX coordinate CSV
   The CSV file should be supplied using;
   `--hox` 

   It must contain "," separated;
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
   Example: [`hox_clusters.annotated.csv`](./hox_clusters.annotated.csv)
```
genome,chrom_file,sequence,hoxgene,scaffold,start,end,strand
Avul.fasta,Avul.chrom,Avul9675_KAN0984125.1,lab,CM147319.1,201577313,201664020,-
Galen.fasta,Galen.chrom,Galen5335_XP_003744805.1,lab,NW_003805462.1,424942,457126,-
Hyas.fasta,Hyas.chrom,Hyas9685_KAH6930615.1,lab,CM023485.1,87184611,87185051,-
Indiac1_clean.fasta,Indi1.chrom,g65748.t1,lab,contig_4489,257705,331956,+
Varr.fasta,Varr.chrom,Varr3255_XP_022670603.1,AbdB,NW_019211460.1,8043840,8047999,+
```
The `genome` column is used to identify the taxa.

# Output
The script generates six files;

```
miRNA_HOX_positions.tsv
miRNA_presence_absence.tsv
Canonical_HOX_miRNAs.tsv
All_miRNA_HOX_distance_matrix.tsv
miRNA_HOX_heatmap.pdf
miRNA_HOX_heatmap.svg
```
### Output files description

 `miRNA_HOX_positions.tsv`

**Contains the positional relationship between every miRNA and HOX cluster.**

Columns include:

- Taxon
- miRNA
- Scaffold
- Position
- Status
- Distance_bp
- HOX_score
- HOX_context

Example:
```
Taxon	miRNA	Scaffold	Position	Status	Distance_bp	HOX_score	HOX_context
species1	mir-10	scaffold_1	102050	INSIDE_HOX	0	10	inside_HoxA1
species1	mir-100	scaffold_1	95000	SAME_SCAFFOLD	5000	6.3001	upstream_of_HoxA1
```
`miRNA_presence_absence.tsv`
**A binary presence/absence matrix with:**

- rows = taxa
- columns = miRNA families
- values = 0 or 1

Example:
```
Taxon	mir-10	mir-100	mir-2
species1	1	1	0
species2	1	0	1
species3	0	1	1
```
`Canonical_HOX_miRNAs.tsv`

**Contains miRNAs whose names match the canonical HOX-associated patterns:**

- Mir-10
- Iab-4
- Mir-2

The results are sorted by taxon and HOX association score.

`All_miRNA_HOX_distance_matrix.tsv`

**Matrix containing the maximum HOX association score for each taxon/miRNA combination.**

Rows represent taxa and columns represent miRNAs.

`miRNA_HOX_heatmap.pdf`

**PDF visualization of the HOX association matrix.**

`miRNA_HOX_heatmap.svg`

**SVG version of the same heatmap, useful for editing in vector graphics software.**

## Mapping categories

- `INSIDE_HOX`
   - The miRNA midpoint falls inside the genomic span of the HOX cluster.
- `SAME_SCAFFOLD`
   - The miRNA is on the same scaffold as the HOX cluster but outside the cluster boundaries.
- `DIFFERENT_SCAFFOLD`
   - The miRNA is not located on a scaffold containing a recognized HOX cluster.

**For miRNAs located on the same scaffold as a HOX cluster, the script determines a more detailed context.**
Which includes;
```
inside_Hox
upstream_of_Hox
downstream_of_Hox
between_HoxA1_and_HoxA2
```
