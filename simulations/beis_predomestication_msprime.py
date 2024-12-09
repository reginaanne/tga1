# msprime coalescent simulation
# generate neutral genetic variation for SLiM forward sim
# demographic model adapted from Beissinger et al 2016
# Outputs ts generated under neutral pre-domestication demography to load into slim that has neutral mutations

##########################################################
# Import libraries
import sys
import tskit
import msprime
import pyslim
import numpy as np

##########################################################
# command line parameters
trees = sys.argv[1]
rec = sys.argv[2]
seed = sys.argv[3]

##########################################################
# demographic model adapted from Beissinger et al 2016
demog_model = msprime.Demography()
demog_model.add_population(initial_size=12278)

# Simulate ancestry for 656 samples (bottleneck population size/starting size for SLiM)
ots = msprime.sim_ancestry(
        samples=656,
        demography=demog_model,
        recombination_rate=rec,
        sequence_length=1000000,
        random_seed=seed)

# Annotate individuals, nodes, etc.
ats = pyslim.annotate(ots, model_type="WF", tick=1, stage="late")

# Add SLiM mutations
mut_model = msprime.SLiMMutationModel(type = 1, slim_generation = 1) # neutral mutations are type 1, negative slim_gen for pre-domestication
ts = msprime.sim_mutations(
							ats,
							rate=3e-7, # mutation rate
							model=mut_model,
							keep=True,
							random_seed=seed
)

# Check that it worked (there should be more mutations)
print(f"The tree sequence now has {ts.num_mutations} mutations, at "
      f"{ts.num_sites} distinct sites.\n"
      f"Mean pairwise nucleotide diversity is {ts.diversity():0.3e}.")

##########################################################
# Final ts prep
# Check metadatatables = ts.tables
tables = ts.tables
ts_metadata = tables.metadata
ts_metadata["SLiM"]["model_type"] = "WF"
tables.metadata = ts_metadata
ts = tables.tree_sequence()

# Change time units to generation
t = ts.dump_tables()
t.time_units = 'generations'
ts2 = t.tree_sequence()

# Outputs the trees
ts2.dump(trees)