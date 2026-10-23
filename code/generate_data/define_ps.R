# Clear environment
rm(list = ls())

library(qiime2R)
library(phyloseq)
library(tidyverse)
source("./code/generate_data/rename_unclassified_seqs.R")

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
  "c_AC" = "Outer",
  "c_ACb" = "Inner"
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

# ------ Save at ASV level ------

saveRDS(ps_filt, file = "./data/phyloseq/ps_ASV.rds")

# ------ Rename Unclassified (NA) sequences and agglomerate ------

# replaces NA genus values with lowest classified taxon
ps_filt2 <- rename_genus(ps_filt) 

# combines taxa with the same genus
ps_genus = tax_glom(ps_filt2, "Genus") 

# relative abundance (%) per sample comprised of unclassified genera
sample_sum <- check_unclassified_sum(ps_genus)

# ------ Save ------

saveRDS(ps_genus, file = "./data/phyloseq/ps_genus.rds")
