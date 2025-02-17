# Calculating and comparing SFS from relate tree sequence with observed VCF
# Not using rate map
# Code adapted from Nate Pope's SINGER Snakemake pipeline

import numpy as np
import tskit
import pandas as pd
import allel
import matplotlib.pyplot as plt

### Tree sequence
# Load in trees
ts = tskit.load("haploB.trees")
mutation_rate = 3.3e-8

# function documentation: https://tskit.dev/tskit/docs/stable/python-api.html#tskit.TreeSequence.allele_frequency_spectrum
afs = ts.allele_frequency_spectrum(mode='branch', span_normalise=True, polarised=True)
afs = afs * mutation_rate

filename = f"ts_branched_folded.csv"
    
header = f"ts_branched_folded"
    
np.savetxt(filename, afs, delimiter=",", header=header, comments="")

### VCF
vcf_file = "haploB_sfs_input.vcf.recode.vcf"
vcf = allel.read_vcf(vcf_file)

# a bunch of filtering steps
genotypes = allel.GenotypeArray(vcf['calldata/GT']) 
positions = vcf['variants/POS']

counts = genotypes.count_alleles(max_allele = 1)

num_samples = 672
min_pos = min(positions)
max_pos = max(positions)
length_chunk = max_pos - min_pos

folded_afs_vcf = allel.sfs_folded(counts, n = 2 * num_samples) / length_chunk

np.savetxt("vcf_folded.csv", folded_afs_vcf, delimiter=",", header="vcf_folded", comments="")

### Plot
freq = np.arange(1, folded_afs_vcf.size)
plt.scatter(freq, afs[1:673], c='firebrick', label='branch-ARG', s=8)
plt.scatter(freq, folded_afs_vcf[1:673], c='black', label='site-VCF', s=8)
plt.xlabel("Minor allele frequency")
plt.ylabel("# of variants / base")
plt.yscale("log")
plt.legend()
plt.tight_layout()
plt.savefig("folded_branchARG_siteVCF.jpg", dpi=300)
plt.clf()