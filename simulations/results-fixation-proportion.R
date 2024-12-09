# Mutation fate output preparation from selection simulations

# Load packages
library(readr)
library(vroom)
library(tidyverse)

# Keep the seeds from
keep_seeds <- read_delim("~/tga1/analysis/slimulate-domestication-sweeps/results/sim_h2h1_300seeds.txt", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)

# Simulation runs to process
sim_runs <- list.dirs(path = "../runs", 
                      full.names = FALSE, 
                      recursive = FALSE)
processed_runs_files <- list.files(path = "../results",
                            pattern = "_fix_generation_all.txt")
processed_runs <- str_remove(processed_runs_files,
                             pattern = "_fix_generation_all.txt")

unprocessed_runs <- sim_runs[!sim_runs %in% processed_runs]

# Mapping for formatting
demographic_model_mapping <- c("beis" = "Beissinger2016",
                               "constant" = "Constant",
                               "msmc" = "Wang2017")
starting_freq_mapping <- c("denovo" = "1/n",
                                "standingvariation1" = "~1%",
                                "standingvariation5" = "~5%")

recombination_rate_mapping_raw <- c("rec6" = "2e-6",
                                   "rec8" = "2e-8")
recombination_rate_mapping_adj <- c("rec6" = "2e-7",
                                    "rec8" = "2e-9")
selection_mapping_raw <- c("sel1" = "0.5",
                           "sel2" = "0.1",
                           "sel3" = "0.05",
                           "sel4" = "0.01")
selection_mapping_adj <- c("sel1" = "0.05",
                           "sel2" = "0.01",
                           "sel3" = "0.005",
                           "sel4" = "0.001")

# Function to process files for each sim run
process_sim_runs <- function(sim_run) {
  directory <- sim_run
  
  # Load in files for lost mutation
  lost_gen_files <- list.files(path = paste0("../runs/", directory),
                               pattern = "_lost_generation.txt",
                               full.names = TRUE)

  data_list_lost_gen <- map(lost_gen_files, ~ {
    vroom(.x, col_names = FALSE) %>%
      select(X1, X6) %>%
      rename(gen = X6, seed = X1) %>%
      mutate(seed = str_remove(seed, ":"))
  })
  
  lost_generation <- bind_rows(data_list_lost_gen)
  
  # Load in files for fixed mutations
  fix_gen_files <- list.files(path = paste0("../runs/", directory),
                              pattern = "_fix_generation.txt",
                              full.names = TRUE)
  
  data_list_fix_gen <- map(fix_gen_files, ~ {
    vroom(.x, col_names = FALSE) %>%
      select(X1, X3, X5, X8) %>%
      rename(seed = X1,
             gen = X3,
             pos = X5,
             start_freq = X8)
  })
  
  fix_generation <- bind_rows(data_list_fix_gen)
  
  # Load in files for segregating mutations
  failed_to_fix_files <- list.files(path = paste0("../runs/", directory),
                                    pattern = "_failedToFix.txt",
                                    full.names = TRUE)
  
  data_list_failed <- map(failed_to_fix_files, ~ {
    vroom(.x, col_names = FALSE) %>%
      select(X1) %>%
      rename(seed = X1) %>%
      mutate(seed = str_remove(seed, ":"))
  })
  
  failed_seeds <- bind_rows(data_list_failed)
  
  # Output compiled information on fixed mutations
  write_delim(fix_generation, 
              file = paste0("../results/", directory, "_fix_generation_all.txt"), 
              delim = "\t")
  
  # Output summary statistics
  # Summary statistics column names preparation
  col_names_summary_stats <- c("min", "first_quartile", "median", "mean", 
                               "third_quartile", "max")
  
  seeds <- c()
  if (nrow(fix_generation) > 0) {
    seeds <- c(seeds, fix_generation$seed)
  }
  
  if (nrow(lost_generation) > 0) {
    seeds <- c(seeds, lost_generation$seed)
  }
  
  if (nrow(failed_seeds) > 0) {
    seeds <- c(seeds, failed_seeds$seed)
  }
  
  # Keep all seeds if there are fewer than 500
  if (length(seeds) < 500) {
    keep_seeds <- seeds
  } else {
    keep_seeds <- sample(seeds, 500)
  }
  
  # Filter and count seeds for each generation
  num_fix <- if (nrow(fix_generation) > 0) {
    fix_generation <- fix_generation %>%
      filter(seed %in% keep_seeds)
    length(fix_generation$seed)
  } else {
    0
  }
  
  num_lost <- if (nrow(lost_generation) > 0) {
    lost_generation <- lost_generation %>%
      filter(seed %in% keep_seeds)
    length(lost_generation$seed)
  } else {
    0
  }
  
  num_seg <- if (nrow(failed_seeds) > 0) {
    failed_seeds <- failed_seeds %>%
      filter(seed %in% keep_seeds)
    length(failed_seeds$seed)
  } else {
    0
  }
  
  total <- sum(num_lost, num_fix, num_seg)
  
  # Summary statistics for fix_generation
  if (nrow(fix_generation) > 0) {
    summary_fix <- summary(fix_generation$gen)
    summary_fix_df_raw <- as.data.frame(as.list(summary_fix),
                                    stringsAsFactors = FALSE,
                                    col.names = col_names_summary_stats) %>%
      mutate(label = directory) %>%
      separate(label, c("demographic_model", "recombination_rate", "starting_freq", "selection"), sep = "_") %>%
      mutate(prop_fix = num_fix/total) %>%
      select(demographic_model, recombination_rate, starting_freq, selection, prop_fix, everything())
    summary_fix_df_adj <- summary_fix_df_raw %>%
      mutate(across(c(min, first_quartile, median, mean, third_quartile, max), ~ . * 10)) %>%
      mutate(demographic_model = recode(demographic_model, !!!demographic_model_mapping),
             recombination_rate = recode(recombination_rate, !!!recombination_rate_mapping_adj),
             starting_freq = recode(starting_freq, !!!starting_freq_mapping),
             selection = recode(selection, !!!selection_mapping_adj))
  }
  
  # Summary statistics for lost_generation
  if (nrow(lost_generation) > 0) {
    summary_lost <- summary(lost_generation$gen)
    summary_lost_df_raw <- as.data.frame(as.list(summary_lost),
                                     stringsAsFactors = FALSE,
                                     col.names = col_names_summary_stats) %>%
      mutate(label = directory) %>%
      separate(label, c("demographic_model", "recombination_rate", "starting_freq", "selection"), sep = "_") %>%
      mutate(prop_lost = num_lost/total) %>%
      select(demographic_model, recombination_rate, starting_freq, selection, prop_lost, everything())
    summary_lost_df_adj <- summary_lost_df_raw %>%
      mutate(across(c(min, first_quartile, median, mean, third_quartile, max), ~ . * 10)) %>%
      mutate(demographic_model = recode(demographic_model, !!!demographic_model_mapping),
             recombination_rate = recode(recombination_rate, !!!recombination_rate_mapping_adj),
             starting_freq = recode(starting_freq, !!!starting_freq_mapping),
             selection = recode(selection, !!!selection_mapping_adj))
  }
  
  # Proportions
  proportion_df_raw <- data.frame(label = directory,
                              mutation_fixed = (num_fix/total),
                              mutation_lost = (num_lost/total),
                              mutation_segregating = (num_seg/total)) %>%
    separate(label, c("demographic_model", "recombination_rate", 
                      "starting_freq", "selection_coefficient"),
             sep = "_")
  
  # Formatting for the names of demographic models, starting freq, selection strength
  proportion_df_adj <- proportion_df_raw %>%
    mutate(demographic_model = recode(demographic_model, !!!demographic_model_mapping),
           recombination_rate = recode(recombination_rate, !!!recombination_rate_mapping_adj),
           starting_freq = recode(starting_freq, !!!starting_freq_mapping),
           selection_coefficient = recode(selection_coefficient, !!!selection_mapping_adj))
    
  # Write outputs
  if (num_fix > 0) {write_delim(summary_fix_df_raw, 
              file = "../results/fix_summary_stats_all_raw.txt", 
              delim = "\t", 
              append = TRUE, 
              col_names = !file.exists("../results/fix_summary_stats_all_raw.txt")) 
  write_delim(summary_fix_df_adj, 
              file = "../results/fix_summary_stats_all_adj.txt", 
              delim = "\t", 
              append = TRUE, 
              col_names = !file.exists("../results/fix_summary_stats_all_adj.txt"))
  }
  
  if (num_lost > 0) {write_delim(summary_lost_df_raw, 
              file = "../results/lost_summary_stats_all_raw.txt", 
              delim = "\t", 
              append = TRUE, 
              col_names = !file.exists("../results/lost_summary_stats_all_raw.txt"))
  write_delim(summary_lost_df_adj, 
              file = "../results/lost_summary_stats_all_adj.txt", 
              delim = "\t", 
              append = TRUE, 
              col_names = !file.exists("../results/lost_summary_stats_all_adj.txt"))
  }
  
  write_delim(proportion_df_raw, 
              file = "../results/proportion_stats_all_raw.txt", 
              delim = "\t", 
              append = TRUE, 
              col_names = !file.exists("../results/proportion_stats_all_raw.txt"))
  write_delim(proportion_df_adj, 
              file = "../results/proportion_stats_all_adj.txt", 
              delim = "\t", 
              append = TRUE, 
              col_names = !file.exists("../results/proportion_stats_all_adj.txt"))
}

# Apply to all
test <- unprocessed_runs[1:5]
map(test, process_sim_runs)
