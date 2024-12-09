#!/bin/bash

test=$1

module load R/4.3.2


printf "\n##########\n starting\n##########\n"
vcftools \
	--vcf /home/reginaf/tga1/analysis/relate-new-vcf-subsets/chrom4_haploB_subset.vcf \
	--recode \
	--stdout \
	--chr 6 \
	--from-bp 84500000 \
	--to-bp 85500000 > $test.vcf

printf "\n##########\n vcf made\n##########\n"

~/src/relate_v1.2.1_x86_64_static/bin/RelateFileFormats \
	--mode ConvertFromVcf \
	--haps $test.haps \
	--sample $test.sams \
	-i ./$test

printf "\n##########\n haps and sams made\n##########\n"

~/src/relate_v1.2.1_x86_64_static/scripts/PrepareInputFiles/PrepareInputFiles.sh \
   	--haps $test.haps \
	--poplabels fullsam.poplabels \
	--sample $test.sams \
	--ancestor diplux.6.fasta \
	-o $test.prepared \
	--mask NP.chr6.fa

printf "\n##########\n other inputs made\n##########\n"

~/src/relate_v1.2.1_x86_64_static/scripts/RelateParallel/RelateParallel.sh \
	-m 3.3E-8 \
	-N 300000 \
	--haps $test.prepared.haps.gz \
	--sample $test.prepared.sample.gz \
	--map ogut.chr6.txt \
	--fb 1 \
	--dist $test.prepared.dist.gz \
	--output $RANDOM \
	--threads 30

printf "\n##########\n anc and mut made\n##########\n"

#~/src/relate_v1.2.1_x86_64_static/scripts/TreeView/TreeViewMutation.sh  \
#	--haps $test.prepared.haps.gz \
#	--sample $test.prepared.sample.gz \
#	--mut $test.mut \
#	--anc $test.anc \
#	--bp_of_interest 46350866 \
#	--years_per_gen 1 -o $test.tree \
#	--poplabels fullsam.poplabels

printf "\n##########\n tree made\n##########\n"
