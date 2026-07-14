suppressPackageStartupMessages({
  # ---- Load Packages ----
  library(tidyverse)
  library(phyloseq)
  # formatting figures
  library(patchwork)
  # colors
  library(RColorBrewer)
  library(MetBrewer) 
})
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE, cache = TRUE)

# load phyloseq object (absolute counts)
ps <- readRDS("./data/phyloseq/ps_genus.rds")

sam_name <- c("inner", "outer")

# Metabolism input file
metab_fname <- "./data/metabolism_midas_mbfr.xlsx"


# Color Palettes (MetBrewer)
size_pal <- "Greek"     # reds
taxa_pal <- "Hiroshige" # orange, blue


# ------ Functions for processing phyloseq object ----

get_metadata <- function(ps){
  # cannot use as.data.frame because class will be phyloseq 
  metadata <- data.frame(sample_data(ps)) %>%
    tibble::rownames_to_column("Sample") %>%
    arrange(size.name)
}

get_taxonomy <- function(ps){
  taxonomy <- data.frame(tax_table(ps)) %>%
    tibble::rownames_to_column("OTU")
}

get_rel <- function(ps) {
  # define relative abundance
  ps_rel <- phyloseq::transform_sample_counts(ps, function(x) x*100/sum(x))
  # combines taxonomy and abundance
  phyloseq::psmelt(ps_rel)
}

get_rel_wide <- function(ps) {
  rel_wide <- get_rel(ps) %>%
    dplyr::select(Genus, Sample, Abundance) %>%  
    pivot_wider(
      names_from = Sample,
      values_from = Abundance
    ) %>%
    column_to_rownames(var = "Genus") %>%
    dplyr::select(all_of(sam_name)) 
}


# input df must contain column named Genus
get_metabolism <- function(df) {
  library(readxl)
  m <- read_excel(metab_fname, sheet = "input") # tibble
  
  metab <- df %>%
    dplyr::select(Genus) %>%
    distinct() %>%
    left_join(., m, by = "Genus") %>%
    column_to_rownames("Genus")
}