# Plotting mutation ages vs. h2/h1

library(readr)
library(vroom)
library(tidyverse)
library(gridExtra)

causative <- 46350866

################################################################################
### H2H1
## Observed
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

## Simulated
simulated_h2h1 <- read_delim("~/tga1/analysis/slimulate-domestication-sweeps/results/sim_h2h1_300seeds.txt", 
                                      delim = "\t", escape_double = FALSE, 
                                      trim_ws = TRUE)

### Allele frequencies (observed)
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

### Mutation ages
## Observed
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

relate_ages_subset <- relate_ages_subset %>%
  bind_cols(freq_subset) %>%
  mutate(frq = if_else(alternative_allele == ALLELE1, ALLELE_FREQ1, ALLELE_FREQ2)) %>%
  mutate(sel_pos = causative,
         sel_age = causative_age,
         label = "empirical") %>%
  rename(pos = pos_of_snp,
         age = age_mid) %>%
  select(pos, age, sel_age, sel_pos, label, frq)

## Simulated
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
  filter(seed %in% simulated_h2h1$seed)

################################################################################
## Add in the h12, h2h1
# Simulated
ages_df_with_h12 <- ages_df_adj %>%
  left_join(simulated_h2h1 %>% select(seed, maize_h12, maize_h2h1), by = "seed") %>%
  select(label, age, sel_age, frq, maize_h12, maize_h2h1) %>%
  rename(h12 = maize_h12,
         h2h1 = maize_h2h1) %>%
  separate(label, into = c("demo", "rec", "freq", "sel"))

# Observed
obs_h12 <- mean(empirical_results_h2h1_subset$maize_h12)
obs_h2h1 <- mean(empirical_results_h2h1_subset$maize_h2h1)

relate_df_with_h12 <- relate_ages_subset %>%
  mutate(h12 = obs_h12,
         h2h1 = obs_h2h1) %>%
  mutate(demo = "observed",
         rec = "observed",
         freq = "observed",
         sel = "observed") %>%
  select(demo, rec, freq, sel, age, sel_age, frq, h12, h2h1)

## Combine observed, simulated into one df
plot_df_h2h1_ages <- rbind(ages_df_with_h12, relate_df_with_h12) %>%
  filter(
    rec %in% c("rec6", "observed"),
    sel %in% c("sel1", "observed"),
    between(frq, 0.05, 0.95)
  ) %>%
  mutate(label = paste(demo, freq, sep = "_"))

################################################################################
plot_df_h2h1_ages %>%
  filter(label != "observed_observed") %>%
  ggplot(aes(x = h2h1, y = log10(age))) +
  geom_density_2d_filled(bins = 20, contour_var = "ndensity") +
  geom_hline(yintercept = log10(15500), linetype = "dashed", size = 1) +
  # scale_fill_manual(values = colorRampPalette(RColorBrewer::brewer.pal(9, "Reds"))(20)) +
  # Observed values
  annotate(
    "point",
    x = filter(plot_df_h2h1_ages, demo == "observed")$h2h1,
    y = log10(filter(plot_df_h2h1_ages, demo == "observed")$age),
    fill = "red",
    color = "white",
    pch = 21,
    size = 4
  ) +
  labs(
    x = "H2H1",
    y = "log10(Age)"
  ) +
  xlim(0, 1) +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_text(size = 16),      
        axis.title = element_text(size = 18),    
        strip.text = element_text(size = 16),     
        plot.title = element_text(size = 20, hjust = 0.5),
        panel.border = element_rect(color = "black", fill = NA, size = 1),
        panel.spacing = unit(15, "pt")) +
  facet_grid(
    rows = vars(demo), 
    cols = vars(freq), 
    labeller = labeller(
      demo = c(
        "observed" = "Observed",
        "beis" = "Beissinger 2016",
        "constant" = "Constant",
        "msmc" = "Wang et al. 2017"
      ),
      freq = c(
        "denovo" = "de novo",
        "standingvariation1" = "~1%",
        "standingvariation5" = "~5%"
      )
    ),
    scales = "fixed"
  )

ggsave("../plots/h2h1_mut_ages_density_rec6.jpg",
       units = "in",
       width = 10,
       height = 10)
