# NOTE: Using all data, not just top n_show
# This is to verify that major metabolic groups are not being dropped

rm(list = ls())
source("./code/01_load_ps.R")
source("./code/02_sum_rel_ab_by_function_mbfr.R")

write2excel <- 0

metab_order <- c("Methanogen", "Acetogen", "Fermentation", "Filamentous")

rel_ab_df <- sum_rel_ab_by_function_mbfr(ps) %>%
  filter(metab_val == "Positive") %>%
  mutate(
    metab = factor(metab, levels = metab_order)
  )

# ------------ Plot ------------------


p <- ggplot(rel_ab_df, aes(x = Sample, y = mean_sum)) + #, fill = metab_val)) +
  geom_col(position = "dodge", width = 0.6, fill = "steelblue") +
  geom_errorbar(
    aes(ymin = mean_sum - sd_sum, ymax = mean_sum + sd_sum),
    width = 0.2,
    position = position_dodge(width = 0.6)
  ) +
  facet_wrap(~metab, scales = "free_y", nrow = 1) +
  labs(
    y = "Relative\nAbundance (%)",
    x = "Region"
  ) +
  # scale_fill_manual(
  #   name = "Functional Group",
  #   values = c("Positive" = "steelblue",
  #              "Positive + Variable" = "lightgray")
  # ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_rect(
      colour = NA # facet label outline
    )
  ) 
  
# Save plot
fname <- "./figures/functional_group_abundance.png"
ggsave(fname, plot = p, width = 6.5, height = 2, dpi = 300)



# ------ Write Data to Excel

if (write2excel == 1) {
  ### Do not exclude data
  # define relative abundance
  rel_wide <- get_rel_wide(ps) %>%
    rownames_to_column(var = "Genus")
  
  new_m <- get_metabolism(rel_wide) %>%
    # true if any metabolic groups in row are defined
    mutate(tf = as.integer(if_any(everything(), ~ !is.na(.x)))) %>%
    rownames_to_column(var = "Genus")
  
  full_df <- left_join(rel_wide, new_m, by = "Genus") %>%
    filter(tf == 1) %>%
    dplyr::select(-tf) %>%
    relocate(where(is.numeric), .after = where(is.character)) %>%
    arrange(Genus)
  
  library(writexl)
  write_xlsx(full_df, path = "./data/functional_rel_ab.xlsx")
}