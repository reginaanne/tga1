library(readr)
library(tidyverse)
library(Hmisc)

# Load in the ratemap
rmap_path <- "ogut_2015_rmap.txt"
rmap <- read_table(rmap_path)

# Estimate recombination rate per bp; cM/Mb to per bp
whole_genome <- rmap %>% 
  # Change in cM divided by change in physical distance; divide by 100 due to centimorgan
  mutate(rec_rate = c(NA,(abs(diff(pos_cM))/diff(pos_bp)/100))) %>%
  filter(!is.na(length_bp) & !is.na(rec_rate)) %>%
  select(chr, pos_bp, length_bp, rec_rate)

# Genome-wide average recombination rate - weighted average by physical distance
weighted_average_rec_rate <- sum(whole_genome$length_bp * whole_genome$rec_rate) / sum(whole_genome$length_bp)
print(weighted_average_rec_rate)

# Weighted quintiles by physical distance
weighted_quintiles <- Hmisc::wtd.quantile(whole_genome$rec_rate, 
                                          weights = whole_genome$length_bp, 
                                          probs = seq(0, 1, 0.2))

mean_rec_rate_by_quintile <- whole_genome %>%
  # Assign each recombination rate to a quintile
  mutate(weighted_quintile = cut(rec_rate, 
                                 breaks = weighted_quintiles, 
                                 include.lowest = TRUE, 
                                 labels = paste("Q", 1:5, sep = ""))) %>%
  # Group recombination rates by quintile and calculate mean recombination rates
  group_by(weighted_quintile) %>%
  summarise(mean_rec_rate = mean(rec_rate, na.rm = TRUE))

# Up recombination rate 10x for simulations
mean_rec_rate_by_quintile %>%
  filter(!is.na(weighted_quintile)) %>%
  rowwise() %>%
  mutate(message = paste0(
    "Mean recombination rate for quintile ", 
    weighted_quintile, 
    " to use in simulation: ", 
    format(10 * mean_rec_rate, scientific = TRUE)
  )) %>%
  pull(message) %>%
  walk(cat, "\n")

# Local recombination rate at tga1
tga1 <- 46350866
tga1_rec_rate <- whole_genome %>%
  mutate(end_bp = pos_bp + length_bp) %>%
  filter(chr == 4,
         pos_bp < tga1 & end_bp > tga1) %>%
  pull(rec_rate)
tga1_rec_rate_sim <- tga1_rec_rate*10
print(paste0("Local tga1 recombination rate to use in simulation: ", tga1_rec_rate_sim))
