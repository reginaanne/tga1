# Formatting MSMC demographic model output for msprime, SLiM
# Model originally from Wang et al. 2017

# Packages
library(readr)
library(dplyr)

# Load Wang et al. model; column names are gen and pop
msmc <- read_delim("msmc_model_original.txt", 
                   delim = "\t", escape_double = FALSE, 
                   trim_ws = TRUE,
                   skip = 1)

msmc <- msmc %>%
  # Round generations and population sizes to whole numbers and downscale 10x
  mutate(gen = round(gen/10),
         pop = round(pop/10)) %>%
  # Fill in for every generation
  complete(gen = seq(min(gen), max(gen), by = 1)) %>%
  fill(pop) %>%
  # Label generation as occurring before or after start of domestication
  mutate(post_domestication = gen <= 1550,
         pre_domestication = gen >= 1550,)

# For SLiM part
dom_gen <- 1550
slim <- msmc %>%
  filter(post_domestication == TRUE) %>%
  # Reverse generations - forward-in-time; first generation = 1
  mutate(forward_gen = abs(gen - 1550) + 1) %>%
  # Retain generations with pop size changes
  group_by(pop) %>%
  slice_min(forward_gen) %>%
  ungroup() %>%
  arrange(forward_gen) 

SLiM_syntax <- paste0(slim$forward_gen, 
                      " early() {\n  p0.setSubpopulationSize(", 
                      slim$pop, ");\n}")


# Output text file to incorporate into SLiM script (manually)
write.table(SLiM_syntax, 
            file = "Wangetal_2017_SLiM_demo_model.txt", 
            eol= "\n",
            row.names = FALSE, 
            col.names=F,
            quote = FALSE)

# For msprime part
msprime <- msmc %>%
  filter(pre_domestication == TRUE) %>%
  # Retain generations with pop size changes
  group_by(pop) %>%
  slice_min(gen) %>%
  ungroup() %>%
  arrange(gen) %>%
  # Generation 0 = start of domestication
  mutate(gen_dom = gen - 1550)

msprime_syntax <- paste("demography.add_population_parameters_change(time=", 
      msprime$gen_dom, ", initial_size=", 
      msprime$pop, ")", sep="")

# Output text file to incorporate into msprime script (manually)
write.table(msprime, 
            file = "Wangetal_2017_msprime_demo_model.txt", 
            eol= "\n",
            row.names = FALSE, 
            col.names=F,
            quote = FALSE)
