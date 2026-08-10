rm(list = ls())
library(phyloseq)
library(tidyverse)
library(vegan)

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
  SampleID = rownames(otu_table),
  biofilm = c(rep("MBfR", 2), rep("AGS", 18))
) 

# ------ NMDS -----


set.seed(2)
nmds <- metaMDS(dist_matrix) # list

nmds_df <- scores(nmds) %>%
  as_tibble(rownames = "SampleID") %>%
  left_join(., metadata, by = "SampleID")

stress_text <- paste0("2D Stress = ", round(nmds$stress,5))
print(stress_text)


# ------ Plot ------

# colors
cols <- c("#8D1C06", "steelblue")

p <- ggplot(nmds_df, aes(NMDS1, NMDS2, color = biofilm)) +
  geom_point() +
  # geom_polygon(alpha = 0.5, aes(fill = biofilm)) +
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  labs(
    x = "Axis 1", 
    y = "Axis 2",
    color = "Biofilm") + 
    # fill = "Biofilm") +
  theme_classic(base_size = 12) 

fname <- "./figures/Figure_3.tif"
ggsave(fname, plot = p, width = 6.5, height = 3, dpi = 300)
