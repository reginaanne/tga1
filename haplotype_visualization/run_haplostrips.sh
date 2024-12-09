#!/bin/bash 

# Haplostrips on polarized VCF

module load R/3.6.3
conda activate haplostrips

VCF=haplostrips.polarized.4.vcf.gz

# 5kb upstream of causative site, plus entire gene 
# ~10Kb centered on causative site

WINDOW=4:46348366-46355118
OUTPUT=strictFilter.polarized.5kb
/home/reginaf/tga1/scripts/haplostrips/haplostrips.py \
    -v $VCF \
    -o $OUTPUT \
    -i $WINDOW \
    -P haplo.meta.txt \
    -p maize,mexicana,parviglumis,teosinte_mixed \
    -r B73_1 \
    -S 1 -t -c 0.05
