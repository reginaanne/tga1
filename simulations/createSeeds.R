library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
demog <- args[1]
rec <- args[2]
freq <- args[3]
sel <- args[4]

filename <- paste(demog, rec, freq, sel, "seeds.txt", sep = "_")

seeds <- runif(n = 200000, min = 0, max = (2^32)-1) %>% 
unique() %>% 
sample(1000) %>% 
format(scientific = FALSE)

seeds <- as.data.frame(format(seeds, scientific = FALSE))

write.table(seeds, filename, sep="\t",
            row.names=FALSE, col.names=FALSE, quote = FALSE)

