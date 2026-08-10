rm(list = ls())
library(phyloseq)
library(tidyverse)
library(vegan)
library(writexl)

# load phyloseq object for all sample sizes
ps_mbfr <- readRDS("./data/phyloseq/ps_ASV_MBfR.rds")
ps_ags  <- readRDS("./data/phyloseq/ps_ASV_AGS.rds")

otu_mbfr <- psmelt(otu_table(ps_mbfr))
otu_ags  <- psmelt(otu_table(ps_ags))

# Combine MBfR and AGS data 
otu_table_long <- bind_rows(otu_mbfr, otu_ags)

otu_table <- pivot_wider(
  otu_table_long, 
  names_from = "OTU", 
  values_from = "Abundance",
  values_fill = 0
  ) %>%
  column_to_rownames("Sample")
  
rarefy_level <- min(sample_sums(ps_mbfr), sample_sums(ps_ags))

set.seed(1)
dist_matrix <- avgdist(otu_table, sample = rarefy_level, iterations = 10, dmethod = "bray")

metadata <- data.frame(
  row.names = rownames(otu_table),
  biofilm = c(rep("MBfR", 2), rep("AGS", 18))
)


# ------ PERMANOVA ------

overall_res <- adonis2(
  dist_matrix ~ biofilm,
  data = metadata,
  permutations = 999
  ) %>%
  rownames_to_column(var = "Data")

# Multivariate homogeneity of groups dispersions
overall_bd <- anova(
  betadisper(dist_matrix, metadata$biofilm)
  )
