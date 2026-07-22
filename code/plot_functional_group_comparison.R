# NOTE: Using all data, not just top n_show
# This is to verify that major metabolic groups are not being dropped

rm(list = ls())
library(ggh4x)
source("./code/01_load_ps.R")
source("./code/02_sum_rel_ab_by_function_mbfr.R")
source("./code/02_sum_rel_ab_by_function_ags.R")

write2excel <- 0

rel_ab_mbfr <- sum_rel_ab_by_function_mbfr(ps) %>%
  filter(metab == "Filamentous") %>%
  rename(exp_category = Sample) %>%
  mutate(biofilm = "MBfR")

# load AGS data
ps_ags <- readRDS("./data/phyloseq/ps_genus_full_AGS.rds")

size <- data.frame(
  ranges = levels(ps_ags@sam_data$size.mm),
  name = levels(ps_ags@sam_data$size.name)
)

rel_ab_ags <- sum_rel_ab_by_function_ags(ps_ags) %>%
  filter(metab == "Filamentous") %>%
  rename(exp_category = size.name) %>%
  mutate(biofilm = "AGS")

rel_ab_df <- bind_rows(rel_ab_mbfr, rel_ab_ags) %>%
  filter(metab_val == "Positive") %>%
  mutate(
    biofilm = factor(biofilm, levels = c("MBfR", "AGS")),
    exp_category = factor(exp_category, levels = c("Inner", "Outer", "Floccular", "S", "M", "L", "XL", "XXL"))
    )

# ------------ Plot ------------------

p <- ggplot(rel_ab_df, aes(x = exp_category, y = mean_sum)) +  # fill = metab_val)) +
  geom_col(position = "dodge", width = 0.6, fill = "steelblue") +
  geom_errorbar(
    aes(ymin = mean_sum - sd_sum, ymax = mean_sum + sd_sum),
    width = 0.2,
    position = position_dodge(width = 0.6)
  ) +
  facet_wrap(~biofilm, scales = "free") +
  force_panelsizes(cols = c(0.4, 1)) +
  labs(
    y = "Relative\nAbundance (%)",
    x = "Biofilm"
  ) +
  # scale_fill_manual(
  #   name = "Functional Group",
  #   values = c("Positive" = "steelblue",
  #              "Positive + Variable" = "lightgray")
  # ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    strip.background = element_rect(
      colour = NA # facet label outline
    )
  ) 
  
# Save plot
fname <- "./figures/functional_group_comparison.png"
ggsave(fname, plot = p, width = 6.5, height = 2.5, dpi = 300)



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