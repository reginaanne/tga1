#!/bin/bash -l
#SBATCH -D /home/reginaf/tga1/analysis/elai-all-chroms/input
#SBATCH -o /home/reginaf/tga1/slurm-log/elai2-stdout-%j.txt
#SBATCH -e /home/reginaf/tga1/slurm-log/elai2-stderr-%j.txt
#SBATCH -J elai
#SBATCH --array=0-9%5
#SBATCH -t 2:00:00
#SBATCH --mem 8G
#SBATCH --partition=low2

# Prepping VCFs to later make bimbams for ELAI

unset CONDA_EXE
module load vcftools
module load bcftools

# Later scripts ID/TAXA interchangeable instead
ID=$1
TAXA=$2
SAMPLES_LIST=${TAXA}_${ID}_elai.txt

VCF_LIST=($(<input_vcf.txt))
VCF=${VCF_LIST[${SLURM_ARRAY_TASK_ID}]}
CHROM=$(echo "$VCF" | sed 's/.*chrom\([0-9]\+\).*/\1/')

vcftools --gzvcf $VCF --keep $SAMPLES_LIST --recode --recode-INFO-all --out ${ID}.${CHROM}.subsample
