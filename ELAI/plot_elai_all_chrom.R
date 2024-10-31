# Plot ELAI results
# Supplemental figures - all chromosomes

library(readr)
library(tidyverse)
library(patchwork)

# Command line - provide sample ID and generations since gene flow specified for ELAI
args <- commandArgs(trailingOnly = TRUE)
sample_spec <- args[1]
gen_spec <- args[2]

# SNP position file paths
pos_files <- list.files(path = "output", pattern = "*snpinfo.txt*", all.files = TRUE)
results_files <- list.files(path = "results", pattern = "_output.txt", all.files = TRUE)

# Allele dosage order
order <- c("maize", "parv", "mex")

# Plot colors
taxa_colors <- c("maize" = "#0A9396", 
                 "mex" = "#AE2012", 
                 "parv" = "#CA6702")

# Plot formatting x-axis in megabases
format_megabases <- function(x) {
  paste0(x / 1e6, "Mb")
}

# Function to prep df
prep_elai <- function(chrom) {
  
  # Filter for the correct SNP position file for this chromosome
  pos_file_keep <- as.data.frame(pos_files) %>%
    separate(pos_files, into = c("sample", "chromosome", "gen", "seed"), sep = "_") %>%
    filter(sample == sample_spec, gen == gen_spec, chromosome == chrom) %>%
    slice(1) %>%
    mutate(filename = paste(sample, chromosome, gen, seed, sep = "_")) %>%
    pull(filename) %>% 
    unique()
  
  # Read the position data
  positions <- vroom(paste0("output/", pos_file_keep),
                     delim = "\t",
                     escape_double = FALSE,
                     trim_ws = TRUE)
  
  # Filter for the correct results file for this chromosome
  results_file_keep <- as.data.frame(results_files) %>%
    separate(results_files, into = c("sample", "chromosome", "gen", "ext"), sep = "_") %>%
    filter(sample == sample_spec, gen == gen_spec, chromosome == chrom) %>%
    mutate(filename = paste(sample, chromosome, gen, ext, sep = "_")) %>%
    pull(filename) %>% 
    unique()
  
  # Read the results data and calculate column means
  df <- read_table(paste0("results/", results_file_keep), col_names = FALSE)
  df_means <- as.data.frame(colMeans(df))
  colnames(df_means) <- "anc"
  
  # Add position and taxa columns
  df_means <- df_means %>%
    filter(row_number() <= n() - 1) %>%  
    mutate(pos = rep(positions$pos, each = 3),
           taxa = rep(order, times = length(positions$pos)),
           chromosome = chrom)  
  
  return(df_means)
  
}

# Function to plot the prepared df
plot_elai <- function(df_means, chrom) {
  p <- df_means %>%
    filter(chromosome == chrom) %>%
    ggplot(aes(x = pos, y = anc, color = taxa)) +
    geom_ribbon(aes(ymin = 0, ymax = anc, fill = taxa),
                show.legend = FALSE) +
    ylim(c(0, 2)) +
    xlab("SNP position") +
    ylab("ancestry allele dosage") +
    ggtitle(sample_spec,
            subtitle = paste0("chromosome ", chrom)) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      panel.grid.major = element_blank(),  
      panel.grid.minor = element_blank(), 
      axis.line = element_line(color = "black"),
      strip.text = element_text()
    ) +
    scale_color_manual(values = taxa_colors) +
    scale_fill_manual(values = taxa_colors) +
    scale_x_continuous(labels = format_megabases,
                       breaks = function(x) pretty(x, n = 5)) +
    scale_y_continuous(expand = c(0, 0)) +
    facet_grid(rows = vars(taxa),
               axes = "all_x")
  
  return(p)
}

# Function to prep and plot for a chromosome
prep_and_plot <- function(chromosome) {
  df_means <- prep_elai(chromosome)
  
  plot_chrom <- plot_elai(df_means, chromosome)
  ggsave(paste0("plots/elai_", sample_spec, "_", chromosome, "_", gens_spec, ".jpg"), 
         plot = plot_chrom, device = "jpg", width = 12, height = 6)
}

# Loop through all chromosomes and generate plots
for (chromosome in 1:10) {
  prep_and_plot(chromosome)
}
