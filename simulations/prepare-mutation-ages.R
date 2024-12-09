# Prepare mutation age results from domestication sweep simulations

################################################################################
# SET UP
################################################################################
# Load packages
library(readr)
library(vroom)
library(tidyverse)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
demo <- args[1] 
rec <- args[2]
freq <- args[3]
sel <- args[4]

# Prep label, paths
label <- paste(c(demo, rec, freq, sel), collapse = "_")

results_path <- "../results"
output_path <- "../runs"

output_file_name <- paste0(label, "_ages_results.txt")

position_file_name <- paste0(label, "_fix_generation_all.txt")

h2h1_file_path <- paste0("../results/", label, "_h2h1_results.txt")

################################################################################
# LOAD IN DATA AND PREP DF
################################################################################
# H2H1 windows
h2h1_windows <- read_delim(h2h1_file_path, 
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE)

# Information on position of selected mutation
positions <- vroom(paste0(results_path, "/", position_file_name),
                   delim = "\t",
                   escape_double = FALSE, 
                   col_names = TRUE) %>%
  distinct()
positions <- positions %>%
  rename(sel_pos = pos)

# Load in mutation age results
temp <- list.files(path = (paste0(output_path, "/", label)),
                   pattern = "_freq_ages.txt",
                   full.names = TRUE)
temp_seeds <- sub(".*/(.*)_freq_ages\\.txt$", "\\1", temp)
keep_files <- paste0(output_path, "/", label, "/", temp_seeds, "_freq_ages.txt")

ages <- lapply(keep_files, function (df) {
  # Extract the file name without the extension to get seed number
  file_name <- tools::file_path_sans_ext(basename(df))
  seed <- stringr::str_extract(file_name, "([^_]+)")
  
  # Read the file using vroom and add seed to the dataframe
  data <- vroom(df,
                delim = ";",
                escape_double = FALSE) %>%
    mutate(seed = seed)
})

ages <- lapply(ages, function(df) {
  df$derived_state <- as.character(df$derived_state) 
  return(df)  
})

ages_df <- bind_rows(ages)

# Get the ages of the selected mutations for each seed
get_sel_ages <- function(seed) {
  selected_pos <- positions %>%
    filter(seed == !!seed) %>%
    pull(sel_pos)
  
  result <- ages_df %>%
    filter(pos == !!selected_pos,
           seed == !!seed)
}

selected_pos_ages_list <- lapply(unique(ages_df$seed),
                                 get_sel_ages)

selected_pos_ages <- do.call(rbind,
                             selected_pos_ages_list) 

# Remove seeds where a mutation stacked on selected site before sweep
remove_stacked_old <- selected_pos_ages %>%
  group_by(seed) %>%
  mutate(num_mut = n()) %>%
  filter(num_mut > 1) %>%
  filter(age > 1551) %>%
  mutate(num_old_mut = n()) %>%
  filter(num_old_mut > 1) %>%
  pull(seed)
remove_stacked_old_seeds <- unique(remove_stacked_old)

# Remove the ages of the stacked mutation that is not the selected mutation
remove_stacked_young <- selected_pos_ages %>%
  group_by(seed) %>%
  mutate(num_mut = n()) %>%
  filter(num_mut > 1) %>%
  filter(age < 1551)

selected_pos_ages <- selected_pos_ages %>%
  filter(!seed %in% remove_stacked_old) %>%
  filter(!(seed %in% remove_stacked_young$seed & age %in% remove_stacked_young$age)) %>%
  rename(sel_pos = pos,
         sel_age = age)

# Add information about that selected mutation age to the df
ages_df2 <- ages_df %>%
  filter(seed %in% selected_pos_ages$seed) %>%
  left_join(selected_pos_ages %>% select(seed, sel_age, sel_pos) %>% distinct(), by = "seed") %>%
  mutate(label = label)

# Keep just those in the h2h1 windows
window_coordinates <- h2h1_windows %>%
  group_by(seed) %>%
  summarise(
    min_start = min(start),
    max_end = max(end)
  ) %>%
  ungroup() %>%
  mutate(seed = as.character(seed))

filtered_ages_df <- ages_df2 %>%
  inner_join(window_coordinates, by = "seed") %>%
  filter(pos >= min_start & pos <= max_end) %>%
  select(-min_start, -max_end, -derived_state)

write.table(filtered_ages_df, file = paste0("../results/", label, "_mut_ages.txt"),
            col.names = TRUE, row.names = FALSE, quote = FALSE)
