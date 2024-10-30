#!/bin/bash -l
#SBATCH -D /home/reginaf/tga1/analysis/elai-all-chroms/input
#SBATCH -o /home/reginaf/tga1/slurm-log/elai2-stdout-%j.txt
#SBATCH -e /home/reginaf/tga1/slurm-log/elai2-stderr-%j.txt
#SBATCH -J elai
#SBATCH --array=1-10
#SBATCH -t 2:00:00
#SBATCH --mem 8G
#SBATCH --partition=high2

# Make bimbam input files to rule ELAI

CHROM=${SLURM_ARRAY_TASK_ID}
TAXA=$1

unset CONDA_EXE
module load plink

TAXA=CIMBL14
plink --vcf "${TAXA}.${CHROM}.subsample.recode.vcf" --recode bimbam --out "${TAXA}_${CHROM}_bimbam.inp"
