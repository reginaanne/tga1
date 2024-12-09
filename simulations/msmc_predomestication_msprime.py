# msprime coalescent simulation
# generate neutral genetic variation for SLiM forward sim
# demographic model adapted from Wang et al. 2017 
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
# demographic model MSMC from Wang et al 2017
demog_model = msprime.Demography()
demog_model.add_population(initial_size=2268)
demog_model.add_population_parameters_change(time=57, initial_size=2624)
demog_model.add_population_parameters_change(time=365, initial_size=3013)
demog_model.add_population_parameters_change(time=731, initial_size=3430)
demog_model.add_population_parameters_change(time=1164, initial_size=3866)
demog_model.add_population_parameters_change(time=1679, initial_size=4310)
demog_model.add_population_parameters_change(time=2291, initial_size=4742)
demog_model.add_population_parameters_change(time=3016, initial_size=5141)
demog_model.add_population_parameters_change(time=3878, initial_size=5489)
demog_model.add_population_parameters_change(time=4901, initial_size=5773)
demog_model.add_population_parameters_change(time=6115, initial_size=5990)
demog_model.add_population_parameters_change(time=7556, initial_size=6147)
demog_model.add_population_parameters_change(time=9268, initial_size=6257)
demog_model.add_population_parameters_change(time=11299, initial_size=6346)
demog_model.add_population_parameters_change(time=13710, initial_size=6513)
demog_model.add_population_parameters_change(time=19972, initial_size=7167)
demog_model.add_population_parameters_change(time=28795, initial_size=9309)
demog_model.add_population_parameters_change(time=41230, initial_size=15308)
demog_model.add_population_parameters_change(time=58754, initial_size=36083)

# simulate ancestry for 2268 samples (starting pop size for SLiM)
ots = msprime.sim_ancestry(
        samples=2268,
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