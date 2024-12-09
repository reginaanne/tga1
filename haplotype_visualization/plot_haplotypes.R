# Plot just maize individuals
# Uses clustering, output from running Haplostrips

# Load packages
library(readr)
library(tidyverse)

################################################################################
# Sample IDs (to subset maize)
keep_samples <- read_delim("haplo.meta.txt", 
                           delim = "\t", escape_double = FALSE, 
                           trim_ws = TRUE)

# Haplostrips
data <- read_table("strictFilter.polarized.5kb.haps")
distances <- read_table("strictFilter.polarized.5kb.distances_tab")

# Prep data
data$POS <- format(data$POS, scientific = FALSE)
data_long <- tidyr::gather(data, key = "Sample", 
                           value = "Value", -`#CHROM`, -POS, -REF, -ALT)
data_long$SampleID <- sub("_.*$", "", data_long$Sample)

sample_order <- rev(distances$haplotype) 
data_long_sorted <- data_long
data_long_sorted$Sample <- factor(data_long_sorted$Sample, 
                                  levels = sample_order)

sort_indices <- order(factor(data_long$Sample, levels = sample_order))

data_long_sorted <- data_long_sorted[sort_indices, ]

data_long_sorted <- data_long_sorted %>%
  left_join(keep_samples %>% select(ID, populations), by = c("SampleID" = "ID"))

data_long_sorted$populations <- as.factor(data_long_sorted$populations)
data_long_sorted$Value <- as.factor(data_long_sorted$Value)

data_long_sorted$ColorLabel <- ifelse(data_long_sorted$Value == "0", "0", 
                                      paste0(data_long_sorted$populations, "_1"))

just_maize <- data_long_sorted %>%
  filter(populations == "maize") %>%
  left_join(distances, by = c("Sample" = "haplotype")) %>%
  mutate(
    ColorLabel = ifelse(
      Value == "0", "0", 
      paste0("distance_", as.character(distance))
    )
  )

just_maize <- just_maize %>%
  mutate(Sample = factor(Sample, levels = unique(Sample)))

color_palette <- c("0" = "#f0f0f0",
                   "distance_0" = "#009688",
                   "distance_1" = "#AD1457",
                   "distance_2" = "#FFC107",
                   "distance_3" = "#FF9800",
                   "distance_4" = "#F44336",
                   "distance_5" = "#448AFF",
                   "distance_6" = "#1565C0",
                   "distane_7" = "#8BC34A",
                   "distance_8" = "#009688",
                   "distance_9" = "#AD1457",
                   "distance_10" = "#FFC107",
                   "distance_11" = "#FF9800",
                   "distance_12" = "#F44336",
                   "distance_13" = "#448AFF",
                   "distance_14" = "#1565C0",
                   "distance_15" = "#8BC34A",
                   "distance_17" = "#009688")


# Plot with updated fill and color scales
p <- ggplot(just_maize, aes(x = POS, y = Sample, fill = ColorLabel, color = ColorLabel)) +
  geom_tile() +
  scale_fill_manual(values = color_palette) +
  scale_color_manual(values = color_palette) +
  labs(x = "SNP", y = "haplotype") +
  theme(
    panel.background = element_rect(fill = "#f0f0f0"),
    legend.position = "none",
    axis.text.y = element_blank(),
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    line = element_blank(),
    axis.ticks.x = element_line(colour = "black"),
    axis.text.x = element_blank(),
    panel.border = element_rect(colour = "black", fill= NA, size = 1.5),
    text = element_text(size = 18)
  ) +
  scale_x_discrete(
    breaks = seq(0, max(data_long_sorted$POS), by = 1000)
  )

# Save
ggsave("maize_haplos.jpg", plot = p, 
       width = 6, height = 6, units = "in", dpi = 300)
