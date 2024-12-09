# msprime/pyslim/tskit processing of ts output from slim
# add on neutral mutations and prep vcf
# main outputs: VCF to get H2H1/H12, mutation ages
# also outputs mutation frequencies

##########################################################
# Import libraries
import sys
import tskit
import msprime
import pyslim
import numpy as np

##########################################################
# Command line parameters
trees = sys.argv[1]
vcf = sys.argv[2]
rec = sys.argv[3]
seed = sys.argv[4]

##########################################################
# Reads in and cleans up tree from SLiM
orig_ts1 = tskit.load(trees)

# Fix old SLiM format
orig_ts = pyslim.update(orig_ts1)

##########################################################
# Check ts and simplify
# Verify everything coalesces - should be 1 root
orig_max_roots = max(t.num_roots for t in orig_ts.trees())
print(f"Maximum number of roots: {orig_max_roots}\n")

## Simplification
rng = np.random.default_rng(seed = int(seed))

# Get individuals that are alive at the end of the SLiM sim
alive_inds = pyslim.individuals_alive_at(orig_ts, 0)

# Keep just 500 individuals that were alive at the end of the simulation
keep_indivs = rng.choice(alive_inds, 500, replace=False)

# Get SLiM nodes for those individuals
keep_nodes = []
for i in keep_indivs:
  keep_nodes.extend(orig_ts.individual(i).nodes)

# Simpify the ts to just have those nodes for the individuals
sts = orig_ts.simplify(keep_nodes, keep_input_roots=True)

# Check that it worked by printing out how many of the individuals there were before and after simplifying
# note that it's fine/expected to have more individuals, because they have nodes that are required to describe the genealogies of the samples
print(f"Before, there were {orig_ts.num_samples} sample nodes (and {orig_ts.num_individuals} individuals)\n"
      f"in the tree sequence, and now there are {sts.num_samples} sample nodes\n"
      f"(and {sts.num_individuals} individuals).")

##########################################################
# Adding neutral mutations
# Mutation ID tracker
next_id = pyslim.next_slim_mutation_id(sts)

# Add the neutral mutations on
mut_model = msprime.SLiMMutationModel(type = 1, slim_generation = 1552) # neutral mutations are type 1, slim_gen to ensure domestication ages positive
ts = msprime.sim_mutations(
           sts,
           rate=3e-7, # mutation rate, rescaled
           model=mut_model,
           keep=True,
           random_seed=seed
)

# Check that it worked (there should be more mutations)
print(f"The tree sequence now has {ts.num_mutations} mutations,\n"
      f"and mean pairwise nucleotide diversity is {ts.diversity():0.3e}.")

##########################################################
# Extract the individuals and output into VCF
# Get the list of "alive" individuals whose nodes are samples
indivlist = []
for i in pyslim.individuals_alive_at(ts, 0):
    ind = ts.individual(i)
    if ts.node(ind.nodes[0]).is_sample():
       indivlist.append(i)
       # if one node is a sample, the other should be also:
       assert ts.node(ind.nodes[1]).is_sample()

# Write VCF
with open(vcf, "w") as vcffile:
    pyslim.convert_alleles(pyslim.generate_nucleotides(ts)).write_vcf(vcffile, individuals=indivlist)

##########################################################
# Output the mutation ages and frequencies
with open(f"{seed}_freq_ages.txt", "w") as f:
    f.write("site.id;pos;derived_state;slim_time;age;count\n")
    for variant in ts.variants():
        genotypes = variant.genotypes
        alleles, counts = np.unique(np.array(variant.alleles)[genotypes], return_counts=True)
        id_site = variant.site.id
        site = ts.site(id_site)
        pos = site.position
        num_mutations = len(site.mutations)
        mutation_data = []
        
        for i in range(num_mutations):
            mutation = site.mutations[i]
            derived_state = mutation.derived_state
            # Initialize a list for slim_times
            slim_times = []
            # Extract all slim_times from the mutation_list
            for mutation_info in mutation.metadata['mutation_list']:
                slim_times.append(mutation_info['slim_time'])
            # Get the maximum slim_time for this mutation
            slim_time = max(slim_times) if slim_times else None
            age = mutation.time
            mutation_info = {
                'derived_state': derived_state,
                'slim_time': slim_time,
                'age': age
            }
            mutation_data.append(mutation_info)
        
        for allele, count in zip(alleles, counts):
            for mutation in mutation_data:
                if str(mutation['derived_state']) == str(allele):
                    f.write(f"{site.id};{pos};{mutation['derived_state']};{mutation['slim_time']};{mutation['age']};{count}\n")
