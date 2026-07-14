rm(list = ls())
library(writexl)
source("./code/01_load_ps.R")

# ------ Option 1: genera above threshold ------

rel_ab_cutoff <- 0.5 # percent

fname_out <- "./data/midas_genera.xlsx"

# define taxa in which at least one sample has abundance > rel_ab_cutoff
taxa_names <- get_rel(ps) %>%
  filter(Abundance > rel_ab_cutoff) %>%
  distinct(Genus) %>%
  pull(Genus)

high_ab_genera <- get_rel_wide(ps) %>%
  rownames_to_column(var = "Genus") %>%
  filter(Genus %in% taxa_names) %>%
  arrange(Genus) %>%
  dplyr::select(Genus)

# ------ Option 2: Top n ------

n_show <- 30

rel_df <- get_rel_wide(ps) %>%
  # Arrange taxa from largest to smallest abundance
  mutate(row_sum = rowSums(.)) %>%
  arrange(desc(row_sum)) %>%
  # Keep the top n_show
  head(., n = n_show) %>%
  dplyr::select(-row_sum) %>%
  rownames_to_column(var = "Genus") %>%
  arrange(Genus) %>%
  dplyr::select(Genus)
 
write_xlsx(rel_df, path = fname_out)
