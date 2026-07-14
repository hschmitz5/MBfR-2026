# Clear environment
rm(list = ls())

library(qiime2R)
library(phyloseq)
library(tidyverse)
source("./code/generate_data/ps_agglom_function.R")

# Import QIIME2 data as phyloseq object
ps_full <- qiime2R::qza_to_phyloseq(
  features = "./data/qiime/table_dada2_ctrdif.qza",
  tree = "./data/qiime/rooted_tree_ctrdif.qza",
  taxonomy = "./data/qiime/taxonomy_ctrdif.qza",
  metadata = "./data/qiime/sample-metadata-ctrdif.tsv"
)

# keep only samples of interest
keep <- sample_names(ps_full) %in% c("c_AC", "c_ACb")
ps <- prune_samples(keep, ps_full)

# Rename samples
region_names <- c(
  "c_AC" = "outer",
  "c_ACb" = "inner"
)

# rename sample names
sample_names(ps) <- region_names[sample_names(ps)]
  
# ------ Filter ------

remove_names <- taxa_names(
  subset_taxa(
    ps,
    Kingdom == "Unassigned" |    
      Order == "Chloroplast" |
      Family == "Mitochondria"
  )
)

keep_taxa <- !(taxa_names(ps) %in% remove_names)

ps_filt <- prune_taxa(keep_taxa, ps)

# ------ Rarefy ------

# define minimum depth to rarefy
rarefy_level <- min(sample_sums(ps_filt))  # lowest number of ASVs per sample

ps_rare <- rarefy_even_depth(
  ps_filt, rarefy_level, rngseed = 1, replace = FALSE, trimOTUs = TRUE, verbose = TRUE
)

# ------ Save at ASV level ------

saveRDS(ps_filt, file = "./data/phyloseq/ps_ASV.rds")
saveRDS(ps_rare, file = "./data/phyloseq/ps_ASV_rarefied.rds")

# ------ Agglomerate, keeping NA values  ------

ps_genus   <- agglom_genus(ps_filt)
ps_species <- agglom_species(ps_filt)

# ------ Save ------

saveRDS(ps_genus, file = "./data/phyloseq/ps_genus.rds")
saveRDS(ps_species, file = "./data/phyloseq/ps_species.rds")
