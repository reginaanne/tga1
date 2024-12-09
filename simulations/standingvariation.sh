#!/bin/bash -l
#SBATCH -o /home/reginaf/tga1/slurm-log/slim-stdout-%j.txt
#SBATCH -e /home/reginaf/tga1/slurm-log/slim-stderr-%j.txt
#SBATCH --array=0-999%100
#SBATCH -J slim
#SBATCH -t 2:00:00
#SBATCH --mem 18G
#SBATCH --partition=low2

# standing variation

# 1) Runs msprime and SLiM to first simulate pre-domestication neutral variation then simulate selective sweep
# 2) Processes tree sequences to add neutral mutations, then outputs mutation ages and allele frequencies
# 3) Processes VCF generated from ts to calculate H2H1/H12 diversity stats

# Parameters
DEMOG=$1 # demographic model: beis, MSMC, or constant
REC=$2 # rescaled recombination rate: 2e-6 or 2e-7
REC_LABEL=$3 # recombination label for files
FREQ_LABEL=$4 # denovo, standingvariation 1 for 1%, or standingvariation5 for 5%
SEL=$5 # selection coefficient
SEL_LABEL=$6 # selection coefficient label for files

LOWER=$7 # lower frequency bound
UPPER=$8 # upper frequency bound

DIR=../runs/${DEMOG}_${REC_LABEL}_${FREQ_LABEL}_${SEL_LABEL}
cd $DIR

# Check if the expected output files exist
if [ -f "${SIMID}_lost_generation.txt" ] || [ -f "${SIMID}_freq_ages.txt" ]; then
    echo "One or both of the files (${SIMID}_lost_generation.txt, ${SIMID}_freq_ages.txt) exist. Exiting."
    exit 1
fi

# Set up array
SEEDS_FILE="${DEMOG}_${REC_LABEL}_${FREQ_LABEL}_${SEL_LABEL}_seeds.txt"
SIMID_LIST=($(<"${SEEDS_FILE}"))
SIMID=${SIMID_LIST[${SLURM_ARRAY_TASK_ID}]} # seed/simulation ID

# Script paths
MSPRIME1=../../sim_scripts/${DEMOG}_predomestication_msprime.py
SLIM=../../sim_scripts/${DEMOG}_standingvariation_domestication_forward.slim
MSPRIME2=../../sim_scripts/end_msprime.py

## SIMULATION AND TREE SEQUENCE PROCESSING
# Predomestication
module load conda/msprime
TREES=${SIMID}_initial.trees
python $MSPRIME1 $TREES $REC $SIMID

# Domestication - forward - sweep
conda activate slim
echo $TREES
slim -d sel=$SEL -d rec=$REC -d lower=$LOWER -d upper=$UPPER -d "trees='$TREES'" -seed ${SIMID} $SLIM
conda deactivate

# Domestication - process tree sequences
python $MSPRIME2 ${SIMID}_slim.trees ${SIMID}.vcf $REC $SIMID

## VCF PROCESSING
module load bcftools
bcftools query -l ${SIMID}.vcf > names_${SIMID}

# Remove non-biallelic SNPs
module load vcftools
vcftools --vcf ${SIMID}.vcf --out filtered_${SIMID} --min-alleles 2 --max-alleles 2 \
--recode --recode-INFO-all

# Add pop column to meta file
awk '{print $0 "\t" "maize"}' names_${SIMID} > meta_${SIMID}

# Calculate the statistics
/home/reginaf/tga1/scripts/lassip/src/lassip --vcf filtered_${SIMID}.recode.vcf \
--hapstats --lassi --winsize 50 --winstep 10 --out ${SIMID} --pop meta_${SIMID} \
--calc-spec --filter-level 2

## CLEAN UP
rm $TREES
rm ${SIMID}.vcf
rm filtered_${SIMID}.recode.vcf
rm names_${SIMID}
rm meta_${SIMID}