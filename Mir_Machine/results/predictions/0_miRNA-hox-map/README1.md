# Hox-miRNA mapping

This script maps miRNA loci to HOX gene clusters across multiple taxa. 
The script reads HOX gene coordinates and miRNA GFF annotations, determines
the genomic relationship between each miRNA and the corresponding HOX cluster, 
calculates an HOX association score, generates presence/absence matrices, extracts 
canonical HOX-associated miRNAs, and produces a heatmap of HOX association scores.

## Requirements
Python 3.8 or newer.

Required Python packages;
pandas
numpy
matplotlib

To install - pip install pandas numpy matplotlib


## Run the script as;

python3 1.hox-miRNA-map.py -m input_of_hox-miRNAs.gffs --hox hox_clusters.annotated.csv

-m input files  --hox hox coordinate CSV file

# Input files

1. miRNA GFF files
   directory should be supplied using;
   -m or --mirna_dir

The directory (gffs_hox-miRNA) should contain files named as:
  filename_miRNA.gff
      Example: Avul_miRNA.gff
               Galen_miRNA.gff
               Ixod_miRNA.gff
               
  The script needs regular 9 column gff/gff3

  Example:
CM147321.1	MirMachine	microRNA	138807572	138807631	50.3	+	.	gene_id=Bantam.PRE;E-value=5.2e-05;sequence_with_30nt=TAACACGACTGGTGGAGGATCAGACAAAACTGGTTTTCACAATGATCATCCAGATGTGTCCGATATCTGAGATCATTGTGAAAGCTGATTTTGTTGTTTCGACAACGAGGGAAGCGGGAC;seed=(p3_seed(GAGATCA*))
CM147319.1	MirMachine	microRNA	208675234	208675291	84.6	+	.	gene_id=Iab-4.PRE;E-value=7.3e-17;sequence_with_30nt=GACGGCGATGGCTGGATCGGCATCTCCTGTTCGTATACTGAGTGTATCCTGAGTGGACAACTTTCCGGTATACCTTCAGTATACGTAACAGGCGACCCGTTTCGAGAGCGACTTGGAT;seed=(p5_seed(CGTATAC*))
CM147320.1	MirMachine	microRNA	73825156	73825227	67.4	+	.	gene_id=Let-7.PRE;E-value=3e-10;sequence_with_30nt=TTAAATCATCGCTTTTGCGTGCGCTCAGTGTGAGGTAGTAGGTTGTATAGTTGAGAACTACATCTGTCTCGGAGGCCTAACTGTACAACTTGCTAACTTACTCTGCGTGGCACTGTTAAACGTCGCCTCTCG;seed=(p5_seed(TGAGGTA,GAGGTAG*))
