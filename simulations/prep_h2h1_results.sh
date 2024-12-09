#!/bin/bash -l
#SBATCH -o /home/reginaf/tga1/slurm-log/SLiM-stdout-%j.txt
#SBATCH -e /home/reginaf/tga1/slurm-log/SLiM-stderr-%j.txt
#SBATCH -J SLiM
#SBATCH -t 1-00:00:00
#SBATCH --mem 8G
#SBATCH --partition=high2

DEMOG=$1
REC_LABEL=$2
FREQ_LABEL=$3
SEL_LABEL=$4

unset CONDA_EXE
module load R
Rscript prep_h2h1_results.R $DEMOG $REC_LABEL $FREQ_LABEL $SEL_LABEL