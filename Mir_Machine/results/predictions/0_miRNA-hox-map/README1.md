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


## Usage

<python3 1.hox-miRNA-map.py -m input_of_hox-miRNAs.gffs --hox hox_clusters.annotated.csv

-m input files  
--hox hox coordinate CSV file

# Input files

## 1. miRNA GFF files
   Directory should be supplied using;
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

## 2. HOX coordinate CSV
   The CSV file should be supplied using;
   --hox 

   It must contain 
   genome	
   chrom_file	
   sequence	
   hoxgene	
   scaffold	
   start	
   end	
   strand
      Example: 
      genome	   chrom_file	sequence	                 hoxgene	scaffold	      start	      end	      strand
      Avul.fasta	Avul.chrom	Avul9675_KAN0984125.1	   lab	   CM147319.1	   201577313	201664020	-
      Galen.fasta	Galen.chrom	Galen5335_XP_003744805.1	lab	   NW_003805462.1	424942	   457126	   -
      Hyas.fasta	Hyas.chrom	Hyas9685_KAH6930615.1	   lab	   CM023485.1	   87184611	   87185051	   -

   
