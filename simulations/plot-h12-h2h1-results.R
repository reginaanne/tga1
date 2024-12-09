### Plot H12/H2H1 results from simulations

################################################################################
# SET UP
################################################################################
## Load packages
library(readr)
library(vroom)
library(tidyverse)
library(ggdist)

## Load in and prepare empirical maize results
# Define the causative mutation
causative <- 46350866

empirical_results_path <- file.path("/home/reginaf/tga1/analysis/h2h1-all-chrom",
                               "h2h1.4.maize.lassip.hap.spectra.gz")

empirical_results <- read_delim(empirical_results_path,
                                delim = "\t",
                                trim_ws = TRUE,
                                skip = 1)

empirical_results_subset <- empirical_results %>%
  # Label to keep track of observed vs. simulated
  mutate(label = "empirical",
         pmut = causative) %>%
  select(start, end, ppos, maize_h12, maize_h2h1, label, pmut) %>%
  # Keep just those windows containing causative mutation
  filter(start < causative & end > causative)

## Load in selective sweep simulation results
# H2H1 results files
sim_h2h1_files <- list.files(path = "../results",
                               pattern = "_h2h1_results.txt",
                               full.names = TRUE)

sim_h2h1_files_labels <- gsub(".*/results/(.*)_h2h1.*", "\\1", sim_h2h1_files)

# Files with information about the mutation position, fixation generation
sim_pos_files <- list.files(path = "../results",
                            pattern = "_fix_generation_all.txt",
                            full.names = TRUE)

sim_pos_files_labels <- gsub(".*/results/(.*)_fix_generation_all.*", "\\1", sim_pos_files)
sim_pos_files <- sim_pos_files[sim_pos_files_labels %in% sim_h2h1_files_labels]

# Load them all in
sim_h2h1 <- lapply(sim_h2h1_files, function(df) {
  vroom(df, delim = "\t",
        escape_double = FALSE,
        trim_ws = TRUE)
})

sim_pos <- lapply(sim_pos_files, function(df){
  vroom(df, delim = "\t",
        escape_double = FALSE,
        trim_ws = TRUE)
})

# Load in the mutation age results for filtering to remove runs that had issues w/ selected pos
temp <- list.files(path = "../results",
                   pattern = "_mut_ages.txt",
                   full.names = TRUE)

ages <- lapply(temp, function(df) {
  data <- vroom(df)
})

ages_df <- do.call(rbind, ages)
keep_seeds <- unique(ages_df$seed)

################################################################################
# FILTERING
################################################################################
crit_gen <- 1000

sim_pos_df <- bind_rows(sim_pos) %>%
  # Remove runs where selected mutation did not pass filtering for mutation ages
  filter(seed %in% keep_seeds) %>%
  # Remove sim runs where the mutation failed to fix by ~5000 years ago (1000 sim gen)
  filter(gen < crit_gen)

sim_h2h1_df <- bind_rows(sim_h2h1) %>%
  filter(seed %in% sim_pos_df$seed)

# Randomly subsample 300 seeds per parameter set
sim_h2h1_df_300 <- sim_h2h1_df %>%
  group_by(label) %>%
  distinct(seed) %>%
  slice_sample(n = 300) %>%
  ungroup() %>%
  inner_join(sim_h2h1_df, by = "seed") %>%
  select(seed, start, end, ppos, maize_h12, maize_h2h1, pmut, label = label.x)

sim_pos_df_300 <- sim_pos_df %>%
  filter(seed %in% sim_h2h1_df_300$seed)

################################################################################
# PREP TO PLOT
################################################################################
## Include just those where the mutation fixes at least 30% of the time
unique_seed_counts <- sim_h2h1_df_300 %>%
  group_by(label) %>%
  summarise(unique_seed_count = n_distinct(seed), .groups = "drop")
labels_keep <- unique_seed_counts %>%
  filter(unique_seed_count == 300) %>%
  pull(label)

sim_h2h1_df_min300 <- sim_h2h1_df_300 %>%
  filter(label %in% labels_keep) %>%
  group_by(seed) %>%
  slice_sample(n = 1) %>%
  ungroup() %>%
  select(seed, pmut, maize_h12, maize_h2h1, label) %>%
  arrange(label) %>%
  separate(label, into = c("demo", "rec", "freq", "sel"), sep = "_")
  
sim_pos_df_min300 <- sim_pos_df_300 %>%
  filter(seed %in% sim_h2h1_df_min300$seed)

################################################################################
# PLOT
################################################################################
# Mappings for labels
demographic_model_mapping <- c("beis" = "Beissinger2016",
                               "constant" = "Constant",
                               "msmc" = "Wang2017")

selection_mapping_adj <- c("sel1" = "0.05",
                           "sel2" = "0.01",
                           "sel3" = "0.005",
                           "sel4" = "0.001")

freq_mapping <- c("denovo" = "1/n",
                  "standingvariation1" = "~1%",
                  "standingvariation5" = "~5%")

# Custom colors for fill (matching freq)
custom_colors <- c("denovo" = "#005F73",   
                   "standingvariation1" = "#CA6702", 
                   "standingvariation5" = "#AE2012") 

# Labels for the x-axis and legend
custom_labels <- c("denovo" = "de novo (1/n)", 
                   "standingvariation1" = "~1%", 
                   "standingvariation5" = "~5%")

# Function to create the plot
create_h2h1_plot <- function(rec, rec_label) {
  
  plot <- sim_h2h1_df_min300 %>%
    # Plot one recombination rate at a time
    filter(rec == !!rec) %>%  
    
    ggplot(aes(x = factor(freq), y = maize_h2h1, fill = factor(freq))) +
    ylim(c(0,1)) +
    
    # Highlight observed values
    annotate("rect", xmin = -Inf, xmax = Inf, 
             ymin = min(empirical_results_subset$maize_h2h1), 
             ymax = max(empirical_results_subset$maize_h2h1),
             fill = "grey80", alpha = 0.3) +
    
    # Half violin plot
    stat_halfeye(
      adjust = 0.5,
      justification = -0.2,
      .width = 0,
      point_color = NA,
      show.legend = FALSE
    ) +
    
    # Boxplot
    geom_boxplot(
      width = 0.12,
      alpha = 0.5,
      outlier.shape = NA
    ) + 
    
    # Dots (rain)
    stat_dots(
      aes(color = factor(freq)),
      side = "left",
      justification = 1.1,
      show.legend = FALSE
    ) +
    
    labs(
      title = "Simulated Domestication Sweep, H2H1 Values",
      subtitle = paste0("Locally Elevated Recombination Rate: ", rec_label),  
      x = NULL, 
      y = "H2H1",
      fill = "Mutation\nStarting Frequency"
    ) +

    # Apply custom colors for fill and dots
    scale_x_discrete(labels = freq_mapping) +
    scale_fill_manual(values = custom_colors, labels = custom_labels) +  
    scale_color_manual(values = custom_colors, labels = custom_labels) + 
    
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      
      strip.text.y.left = element_text(angle = 0, hjust = 1, size = 12),  
      strip.text.x = element_text(size = 12), 
      
      plot.title = element_text(size = 16),
      plot.subtitle = element_text(size = 12),
      
      axis.text.x = element_text(size = 10),
      axis.text.y = element_blank(),
      axis.line = element_line(color = "black"),

      strip.placement = "outside"
    ) +
    
    facet_grid(rows = vars(sel), cols = vars(demo),
               labeller = as_labeller(c(demographic_model_mapping, selection_mapping_adj)),
               switch = "y") +  
    coord_flip()
  
  # Return the plot object
  return(plot)
}

# Make and save the plots
plot_rec6 <- create_h2h1_plot(
  rec = "rec6",          
  rec_label = "1e-7"    
)

ggsave("supp_sim_h2h1_recLDhelmet.jpg",
       path = "../plots",
       plot = plot_rec6,
       device = "jpg",
       height = 8, width = 12, units = "in")

plot_rec8 <- create_h2h1_plot(
  rec = "rec8",          
  rec_label = "1e-9"    
)

ggsave("supp_sim_h2h1_recOgut.jpg",
       path = "../plots",
       plot = plot_rec8,
       device = "jpg",
       height = 8, width = 12, units = "in")

write.table(sim_h2h1_df_min300, 
            file = "../results/sim_h2h1_300seeds.txt", 
            sep = "\t", 
            quote = FALSE,
            row.names = FALSE, 
            col.names = TRUE)  

