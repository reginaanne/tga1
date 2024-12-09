#!/bin/bash

# (1) Prepares a diploperennis "consensus" sequence to polarize the VCF with
# (2) Polarizes VCF
# (3) Prepares VCF for running haplostrips

module load vcftools
module load bcftools
module load tabix

ORIGVCF=/home/reginaf/tga1/data_filtered/chrom4_tga1_geno_phased.vcf.gz
DIPLOSAMPLES=diploperennis_samples
ZEAMAYSSAMPLES=mays_samples

###
# Subset VCF to keep just diploperennis individuals
vcftools --gzvcf $ORIGVCF --keep $DIPLOSAMPLES --out diplo --recode --recode-INFO-all

# Use custom perl script to prepare consensus sequence for diploperennis
perl dipcons.pl diplo.recode.vcf 0.8
bgzip diplo.recode.vcf.consensus.vcf
tabix diplo.recode.vcf.consensus.vcf.gz

# Subset VCF to keep maize, mex, parv, and teosinte mixed individuals
vcftools --gzvcf $ORIGVCF --keep $ZEAMAYSSAMPLES --out not_diplo --recode --recode-INFO-all &
bgzip not_diplo.recode.vcf
tabix not_diplo.recode.vcf.gz

# Merge the diploperennis consensus with the other samples
bcftools merge not_diplo.recode.vcf.gz diplo.recode.vcf.consensus.vcf.gz > merged.vcf

# Remove invariant sites
bcftools view -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES' merged.vcf > merged_filtered.vcf

bcftools query -l merged_filtered.vcf > samples_list

###
# Polarizes VCF by the diploperennis consensus
python polarizeVCFbyOutgroup.py -vcf merged_filtered.vcf -out chr4_polarized.vcf -ind 673 -add

###
# Remove diploperennis consensus from VCF
vcftools --vcf chr4_polarized.vcf --keep $ZEAMAYSSAMPLES --out remove_diplo --recode --recode-INFO-all

# Remove invariant sites
bcftools view -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES' remove_diplo.recode.vcf > haplostrips.polarized.4.vcf
bgzip haplostrips.polarized.4.vcf
