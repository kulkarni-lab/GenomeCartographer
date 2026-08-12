# Plot hox-miRNA #

This script generates comparative genomic maps showing the positions of HOX genes and selected HOX-associated miRNAs across multiple taxa.

The script is designed as a visualization tool following the [`1.hox-miRNA-map.py`](./1.hox-miRNA-map.py) script.

## Requirements
Python 3.8 or newer is recommended.

Required Python packages;
`pandas` `matplotlib` `numpy`

To install;
```
pip install pandas matplotlib numpy 
```


## Usage
```
python3 2.plot-hox-miRNA.py --hox hox_clusters.annotated.csv --mir Canonical_HOX-miRNAs.tsv --mirnas Mir-10.PRE Iab-4.PRE
```
`--hox` path to hox coordinate CSV file

`--mir` path to canonical HOX miRNA TSV file

`--mirnas` one or more miRNAs to visualize (Default: Mir-10.PRE, Iab-4.PRE)

# Input files

## 1. HOX coordinate file
   Hox coordinate CSV file should be supplied using;
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
## 2. Canonical HOX miRNA table
  The miRNA table is supplied using;
  
  `--mir`  **Eg:** **--mir Canonical_HOX_miRNAs.tsv**

  This file is one of the output generated from previous HOX-miRNA mapping script [`README1_Hox-miRNA-map.md`](./README1_Hox-miRNA-map.md)

It should contain columns including;
```
Taxon
miRNA
Scaffold
Position
Status
Distance_bp
HOX_score
HOX_context
```

Example: 
```
Taxon       miRNA          Scaffold     Position     Status       Distance_bp   HOX_score   HOX_context
species1   Mir-10.PRE     scaffold_1     102000    INSIDE_HOX        0             10       inside_HoxA1
species1   Iab-4.PRE      scaffold_1     95000     SAME_SCAFFOLD    5000          6.30      upstream_of_HoxA1
species2   Mir-10.PRE     scaffold_7     53000     INSIDE_HOX        0             10       inside_HoxA2
```

## 3. Individual miRNAs to plot
  The argument accepts multiple miRNA names separated by spaces. The default is:

  `Mir-10.PRE` `Iab-4.PRE`

  **Example:** `--mirnas Mir-10.PRE Iab-4.PRE`

# Output 

The output will be the genomic map figure in .PDF and .SVG format.
```
Written:
  Panel1A_Hox_miRNA_map.pdf
  Panel1A_Hox_miRNA_map.svg
```

### Filtering and visualization

The script first selects only the requested miRNAs, it retains miRNAS with either `INSIDE_HOX` or `SAME_SCAFFOLD`. 

miRNAs classified as: `DIFFERENT_SCAFFOLD` are excluded from the visualization.

Each `HOX gene` is represented by _rectangle_ The rectangle is positioned according to the gene's genomic coordinates relative to the minimum and maximum coordinates represented on that scaffold. 

`HOX gene` names are displayed above the _genomic blocks_. `Gene labels` are rotated by _90 degrees_ to reduce overlap.

Selected miRNAs are represented using different marker shapes;
```
Mir-10.PRE	            Circle
Iab-4.PRE	              Square
Other selected miRNAs	  Triangle
```

[!NOTE] 
- The script only plots `INSIDE_HOX`, `SAME_SCAFFOLD`, IT DOES NOT PLOT: `DIFFERENT_SCAFFOLD`.
- The script uses fixed font sizes and labels thus, there are chances of labels to **OVERLAP**
- The figure is a comparative genomic visualization rather than a scale-preserving whole-genome plot.
- Absolute distances should be obtained from `miRNA_HOX_positions.tsv`.

## Recommended workflow
This script is intended to be used after the [`1.hox-miRNA-map.py`](./1.hox-miRNA-map.py) script.

```
Step 1: Map miRNAs to HOX regions
python3 1.hox-miRNA-map.py \
    -m gffs_hox-miRNA \
    --hox hox_clusters.annotated.csv
```

This generates:

`Canonical_HOX_miRNAs.tsv`

among other outputs.
```
Step 2: Generate the genomic map

Use the canonical miRNA table:

python3 3.plot.py \
    --hox hox_clusters.annotated.csv \
    --mir Canonical_HOX_miRNAs.tsv \
    --mirnas Mir-10.PRE Iab-4.PRE Mir-279.PRE Mir-2.PRE
```
