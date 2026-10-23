# If a taxon is unclassified (NA), then it is renamed based on the lowest classified taxon.
# All unclassified taxons sharing that lowest classified taxon are then combined into the same entry.
rename_genus <- function(ps) {
 
  # If Genus is NA, then replace with lowest classified taxon 
  taxonomy <- data.frame(ps@tax_table) %>%
    rownames_to_column("OTU") %>%
    mutate(
      Genus = case_when(
        !is.na(Genus)  ~ Genus,
        !is.na(Family) ~ ifelse(startsWith(Family, "midas"), paste0("Unk_", Family), paste0("Unk_f_", Family)),
        !is.na(Order)  ~ ifelse(startsWith(Order, "midas"), paste0("Unk_", Order), paste0("Unk_o_", Order)),
        !is.na(Class)  ~ ifelse(startsWith(Class, "midas"), paste0("Unk_", Class), paste0("Unk_c_", Class)),
        !is.na(Phylum) ~ ifelse(startsWith(Phylum, "midas"), paste0("Unk_", Phylum), paste0("Unk_p_", Phylum)),
        is.na(Phylum)  ~ paste0("Unknown_Phylum"),
        .default = Genus
      )
    ) %>%
    column_to_rownames("OTU") %>%
    as.matrix()
  
  ps_new <- phyloseq(
    otu_table(ps),
    sample_data(ps),
    tax_table(taxonomy), # new names
    phy_tree(ps)            
  )
  
}

# determines relative abundance of unclassified sequences per sample
check_unclassified_sum <- function(ps) {
  sample_sum <- transform_sample_counts(ps_genus, function(x) x*100/sum(x)) %>%
    psmelt(.) %>% # combines taxonomy and abundance
    dplyr::select(Genus, Sample, Abundance) %>%  
    pivot_wider(
      names_from = Sample,
      values_from = Abundance
    ) %>%
    filter(str_detect(Genus, "Unk")) %>%
    column_to_rownames(var = "Genus") %>%
    colSums()
}
