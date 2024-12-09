# Outputs H2H1/H12 values for windows containing selection mutation
# Removes simulations where the mutation did not fix before 5000 years ago

## Set up
# Load packages
library(readr)
library(vroom)
library(tidyverse)
library(dplyr)

# Command line arguments
args <- commandArgs(trailingOnly = TRUE)
demo <- args[1]
rec <- args[2]
startingFreq <- args[3]
selCoeff <- args[4]

# Prep label, paths
label <- paste(c(demo, rec, startingFreq, selCoeff), collapse = "_")

results_path <- "../results"
output_path <- "../runs"

position_file_name <- paste0(label, "_fix_generation_all.txt")

output_file_name <- paste0(label, "_h2h1_results.txt")

## Load in and prepare data
# List files with the H12, H2H1 statistics
temp <- list.files(path = (paste0(output_path, "/", label)),
                   pattern = ".spectra.gz",
                   full.names = TRUE)

# Load in all of the files
spectra <- lapply(temp, function(df) {
  # extract the file name without the extension to get seed number
  file_name <- gsub("\\.maize.lassip\\.hap\\.spectra\\.gz$", "", basename(df))
  
  # read the data from the file using vroom
  data <- vroom(df, delim = "\t", skip = 1, escape_double = FALSE, trim_ws = TRUE)
  
  # add seed number
  data$seed <- file_name
  
  return(data)
})

# Load in mutation position data
positions <- vroom(paste0(results_path, "/", position_file_name),
                   delim = "\t",
                   escape_double = FALSE, 
                   col_names = TRUE) %>%
  distinct()

## Prep results
# Pull out the relevant columns
results <- lapply(spectra, function(df) {
  subset(df,
         select = c("seed", "start", "end",
                    "ppos", "maize_h12", "maize_h2h1"))
})

# Function to pull out just the windows containing the mutation
subset_mut <- function(df) {
  seed_value <- unique(df$seed)
  
  # Extract mutations for the given seed
  mutation <- positions %>%
    filter(seed == seed_value) %>%
    pull(pos)
  
  # Check if mutation is not empty
  if (length(mutation) > 0) {
    # Modify the dataframe
    df <- df %>%
      mutate(pmut = mutation) %>%
      filter(start <= pmut, end >= pmut)
  }
  
  return(df)
}

results_subset <- lapply(results,
                         subset_mut)
results_df <- bind_rows(results_subset)

# Add label
results_df <- results_df %>%
  mutate(label = label)

# Save output
options(scipen = 999, digits = 10)

output_file <- file.path(results_path,
                         output_file_name)
write.table(results_df,
            file = output_file,
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)
