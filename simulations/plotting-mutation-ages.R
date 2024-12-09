### Plotting mutation age results from simulated domestication sweeps

################################################################################
# SET UP
################################################################################
library(readr)
library(vroom)
library(tidyverse)

causative <- 46350866

## OBSERVED/EMPIRICAL - ESTIMATES FROM RELATE
# Load in H2H1 for the window ranges
empirical_results_path <- file.path("/home/reginaf/tga1/analysis/h2h1-all-chrom",
                                    "h2h1.4.maize.lassip.hap.spectra.gz")
empirical_results_h2h1 <- read_delim(empirical_results_path,
                                delim = "\t",
                                trim_ws = TRUE,
                                skip = 1)
empirical_results_h2h1_subset <- empirical_results_h2h1 %>%
  filter(start < causative & end > causative)
start_tga1_windows <- min(empirical_results_h2h1_subset$start)
end_tga1_windows <- min(empirical_results_h2h1_subset$end)

# Load in empirical/observed mutation age estimates
relate_ages_path <- "/home/jri/projects/ibd/relate/data/haploB_tga1/haploB.mut"
relate_ages <- read_delim(relate_ages_path,
                          delim = ";", escape_double = FALSE, trim_ws = TRUE) %>%
  # Remove those that did not map correctly
  filter(is_not_mapping == 0) %>%
  # Get midpoint to use for the ages
  mutate(age_mid = (age_end + age_begin)/2)

relate_ages_subset <- relate_ages %>%
  filter(pos_of_snp > start_tga1_windows & pos_of_snp < end_tga1_windows) %>%
  separate(`ancestral_allele/alternative_allele`, into = c("ancestral_allele", "alternative_allele"), sep = "/")

# Subset to grab mutation ages within H2H1 windows for empirical relate estimates
causative_age <- relate_ages %>%
  filter(pos_of_snp == causative) %>%
  pull(age_mid)

# Allele frequencies for empirical maize
freq <- vroom("~/tga1/analysis/relate-new-vcfs/maize_allele_frq.frq", 
                   delim = "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  separate(`{ALLELE:FREQ}`, into = c("ALLELE", "FREQ"), sep = "\t") %>%
  separate(ALLELE, into = c("ALLELE1", "ALLELE_FREQ1"), sep = ":") %>%
  separate(FREQ, into = c("ALLELE2", "ALLELE_FREQ2"), sep = ":") %>%
  mutate(ALLELE_FREQ1 = as.numeric(ALLELE_FREQ1),
         ALLELE_FREQ2 = as.numeric(ALLELE_FREQ2))
freq_subset <- freq %>%
  filter(POS %in% relate_ages_subset$pos_of_snp) %>%
  select(ALLELE1, ALLELE_FREQ1, ALLELE2, ALLELE_FREQ2)

relate_ages_subset <- relate_ages_subset %>%
  bind_cols(freq_subset) %>%
  mutate(frq = if_else(alternative_allele == ALLELE1, ALLELE_FREQ1, ALLELE_FREQ2)) %>%
  mutate(sel_pos = causative,
         sel_age = causative_age,
         label = "empirical") %>%
  rename(pos = pos_of_snp,
         age = age_mid) %>%
  select(pos, age, sel_age, sel_pos, label, frq)

## SIMULATED
# Matching seeds from the H2H1/H12 plots
seeds <- read_delim("~/tga1/analysis/slimulate-domestication-sweeps/results/sim_h2h1_300seeds.txt", 
                    delim = "\t", escape_double = FALSE, 
                    trim_ws = TRUE)

# Load in the simulated results
temp <- list.files(path = "../results",
                   pattern = "_mut_ages.txt",
                   full.names = TRUE)
ages <- lapply(temp, function(df) {
  data <- vroom(df)
})

ages_df <- do.call(rbind, ages)
ages_df2 <- ages_df %>%
  mutate(frq = as.numeric(count/1000)) %>%
  select(pos, age, sel_age, sel_pos, label, frq, seed)

# Remove mutation where an old one was fixed and a new one arose on top of it
ages_remove <- ages_df2 %>%
  group_by(seed, pos) %>%
  mutate(num_mut = n(),
         total_frq = sum(frq)) %>%
  filter(num_mut > 1, total_frq == 1) %>%
  filter(age == max(age)) %>%
  ungroup()

ages_df3 <- ages_df2 %>%
  anti_join(ages_remove, by = c("seed", "pos"))

# Adjust age for the downscaling needed to run mutations
ages_df_adj <- ages_df3 %>%
  mutate(age = age * 10,
         sel_age = sel_age * 10) %>%
  filter(seed %in% seeds$seed) %>%
  select(-seed)

# Combine simulated results with the empirical ones
data <- rbind(ages_df_adj, relate_ages_subset) %>%
  mutate(mut_older = sel_age > age) %>%
  mutate(label_copy = label) %>%  
  separate(label_copy, into = c("demo", "rec", "freq", "sel"), fill = "right") %>%
  mutate(across(c(rec, freq, sel), ~ coalesce(., demo))) 

################################################################################
# PLOTTING
################################################################################
# Prep to include observed mutation age estimates in each
empirical_data <- data %>%
  filter(demo == "empirical") %>%
  select(-rec, -freq)

len_data <- length(empirical_data$pos)

recs <- c("rec6", "rec8")
freqs <- c("denovo", "standingvariation1", "standingvariation5")
rec_col <- rep(recs, each = (3 * len_data))
freqs_col <- rep(rep(freqs, each = (len_data)), times = 2)

empirical_data_rep <- bind_rows(replicate(6, empirical_data, simplify = FALSE))

empirical_data_rep$rec <- rec_col
empirical_data_rep$freq <- freqs_col

data2 <- data %>%
  filter(label != "empirical") %>%
  rbind(empirical_data_rep)

data_plot <- data2 %>%
  # Remove intermediate frequency mutations
  filter(between(frq, 0.05, 0.95)) %>%
  
  # Remove the age of the selected mutation itself
  filter(pos != sel_pos) %>%
  
  # Thin down to reduce overplotting
  group_by(label) %>%
  slice_sample(n = 500) %>%  
  ungroup()

# Prep labels
demographic_model_mapping <- c("Observed",
                               "beis" = "Beissinger2016",
                               "constant" = "Constant",
                               "msmc" = "Wang2017")

selection_mapping_adj <- c("sel1" = "0.05",
                           "sel2" = "0.01",
                           "sel3" = "0.005",
                           "sel4" = "0.001")

freq_mapping <- c("denovo" = "1/n",
                  "standingvariation1" = "~1%",
                  "standingvariation5" = "~5%")

# Define custom colors for fill (matching demo)
custom_colors <- c("Observed" = "#001219",
                   "Constant" = "#0A9396",   
                   "Beissinger2016" = "#CA6702", 
                   "Wang2017" = "#AE2012") 

# Define labels for the x-axis and legend
custom_labels <- c("observed" = "observed",
                   "denovo" = "de novo (1/n)", 
                   "standingvariation1" = "~1%", 
                   "standingvariation5" = "~5%")

# Define labels for recombination rate
rec_mapping <- c("rec6" = "1e-7",
                 "rec8" = "1e-9")

# Create the plot
ages_plot <- data_plot %>%
  # Subset for this particular plot
  filter(sel %in% c("sel1", "empirical")) %>%
  
  # Map demographic model levels
  mutate(demo = factor(demo, 
                       levels = c("empirical", "beis", "constant", "msmc"),
                       labels = c("Observed", "Beissinger2016", "Constant", "Wang2017"))) %>%
  
  # Create the plot
  ggplot(aes(x = demo, y = log10(age))) +
  
  # Plot
  geom_jitter(aes(color = demo), alpha = 0.2, height = 0) +
  geom_boxplot(aes(fill = demo), outlier.shape = NA) +
  
  # Line at time selection started
  geom_hline(yintercept = log10(15500), color = "grey40", linetype = "dashed") +
  
  # Labels
  labs(
    title = "Mutation Ages Within Simulated Domestication Sweep (s = 0.05)",
    subtitle = "Mutation Starting Frequency",
    y = "log10 generation in the past",
    x = "Demographic Model",
    color = "Demographic Model",
    fill = "Demographic Model"
  ) +
  
  # Theme
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    
    strip.text.y = element_text(size = 14),  
    strip.text.x = element_text(size = 14),
    
    axis.text.x = element_blank(),
    axis.text.y = element_text(color = "black", size = 14),
    axis.line = element_line(color = "black"),
    plot.title = element_text(size = 18),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    
    legend.text = element_text(size = 12),  
    legend.title = element_text(size = 14)
  ) +
  
  # Set custom colors
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  
  # Facet with custom labels
  facet_grid(cols = vars(freq),
             rows = vars(rec),
             labeller = as_labeller(c(custom_labels, rec_mapping)))

# Plot the result
ages_plot

ggsave("sel1_ages_all.jpg",
       plot = ages_plot,
       path = "../plots",
       width = 12,
       height = 8,
       units = "in")

### For supplemental
# Make separate plots for additional selection coefficients
ages_plot_sel2 <- data_plot %>%
  # Subset for this particular plot
  filter(sel %in% c("sel2", "empirical"),
         freq != "denovo") %>%
  
  # Map demographic model levels
  mutate(demo = factor(demo, 
                       levels = c("empirical", "beis", "constant", "msmc"),
                       labels = c("Observed", "Beissinger2016", "Constant", "Wang2017"))) %>%
  
  # Create the plot
  ggplot(aes(x = demo, y = log10(age))) +
  
  # Plot
  geom_jitter(aes(color = demo), alpha = 0.2, height = 0) +
  geom_boxplot(aes(fill = demo), outlier.shape = NA) +
  
  # Line at time selection started
  geom_hline(yintercept = log10(15500), color = "grey40", linetype = "dashed") +
  
  # Labels
  labs(
    title = "Mutation Ages Within Simulated Domestication Sweep (s = 0.01)",
    subtitle = "Mutation Starting Frequency",
    y = "log10 generation in the past",
    x = "Demographic Model",
    color = "Demographic Model",
    fill = "Demographic Model"
  ) +
  
  # Theme
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    
    strip.text.y = element_text(size = 14),  
    strip.text.x = element_text(size = 14),
    
    axis.text.x = element_blank(),
    axis.text.y = element_text(color = "black", size = 14),
    axis.line = element_line(color = "black"),
    plot.title = element_text(size = 18),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    
    legend.text = element_text(size = 12),  
    legend.title = element_text(size = 14)
  ) +
  
  # Set custom colors
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  
  # Facet with custom labels
  facet_grid(cols = vars(freq),
             rows = vars(rec),
             labeller = as_labeller(c(custom_labels, rec_mapping)))

ages_plot_sel2

ggsave("sel2_ages_all.jpg",
       plot = ages_plot_sel2,
       path = "../plots",
       width = 12,
       height = 8,
       units = "in")

ages_plot_sel3 <- data_plot %>%
  # Subset for this particular plot
  filter(sel %in% c("sel3", "empirical"),
         freq != "denovo") %>%
  
  # Map demographic model levels
  mutate(demo = factor(demo, 
                       levels = c("empirical", "beis", "constant", "msmc"),
                       labels = c("Observed", "Beissinger2016", "Constant", "Wang2017"))) %>%
  
  # Create the plot
  ggplot(aes(x = demo, y = log10(age))) +
  
  # Plot
  geom_jitter(aes(color = demo), alpha = 0.2, height = 0) +
  geom_boxplot(aes(fill = demo), outlier.shape = NA) +
  
  # Line at time selection started
  geom_hline(yintercept = log10(15500), color = "grey40", linetype = "dashed") +
  
  # Labels
  labs(
    title = "Mutation Ages Within Simulated Domestication Sweep (s = 0.005)",
    subtitle = "Mutation Starting Frequency",
    y = "log10 generation in the past",
    x = "Demographic Model",
    color = "Demographic Model",
    fill = "Demographic Model"
  ) +
  
  # Theme
  theme_minimal() +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    
    strip.text.y = element_text(size = 14),  
    strip.text.x = element_text(size = 14),
    
    axis.text.x = element_blank(),
    axis.text.y = element_text(color = "black", size = 14),
    axis.line = element_line(color = "black"),
    plot.title = element_text(size = 18),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    axis.title.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    
    legend.text = element_text(size = 12),  
    legend.title = element_text(size = 14)
  ) +
  
  # Set custom colors
  scale_color_manual(values = custom_colors) +
  scale_fill_manual(values = custom_colors) +
  
  # Facet with custom labels
  facet_grid(cols = vars(freq),
             rows = vars(rec),
             labeller = as_labeller(c(custom_labels, rec_mapping)))

ages_plot_sel3

ggsave("sel3_ages_all.jpg",
       plot = ages_plot_sel3,
       path = "../plots",
       width = 12,
       height = 8,
       units = "in")
