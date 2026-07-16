source activate
conda activate /pfs/home/siddharth/miniconda3/envs/singularity

ASSEMBLY1=/path/genomes/ACARI_Argas_vulgaris.fna
singularity exec  --bind /pfs/  --cleanenv /path/dfam-tetools-latest.sif RepeatModeler -database ACARI_Argas_vulgaris  -threads 30 -LTRStruct
singularity exec  --bind /pfs/  --cleanenv /path/dfam-tetools-latest.sif RepeatMasker -lib ACARI_Argas_vulgaris-families.fa $ASSEMBLY1  -pa 27 -xsmall