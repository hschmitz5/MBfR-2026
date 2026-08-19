rm(list = ls())
library(tidyverse)
library(ggh4x)

import_data <- function(PN_fname, PS_fname, biofilm_type, group_name) {
  PN <- readRDS(PN_fname) %>%
    select(extract, group_name, replicate, PN = C_TSS) 
  
  PS <- readRDS(PS_fname) %>%
    select(extract, group_name, replicate, PS = C_TSS)
  
  df_joined <- left_join(PN, PS, by = c("extract", group_name, "replicate")) %>%
    mutate(
      biofilm = biofilm_type
    ) %>%
    select(biofilm, extract, exp_group = group_name, replicate, PN, PS)
}

ags_data <- import_data("./data/EPS/PN_conc_ags.rds", "./data/EPS/PS_conc_ags.rds", "AGS", "size") 

mbfr_data <- import_data("./data/EPS/PN_conc_mbfr.rds", "./data/EPS/PS_conc_mbfr.rds", "Methanogenic", "region") %>%
  mutate(replicate = as.character(replicate))

# Combine both data sets
eps <- bind_rows(ags_data, mbfr_data)

summary_wide <- eps %>%
  group_by(biofilm, extract, exp_group) %>%
  summarize(
    # protein
    PN_avg = mean(PN, na.rm = TRUE), PN_sd = sd(PN, na.rm = TRUE),
    # polysaccharide
    PS_avg = mean(PS, na.rm = TRUE), PS_sd = sd(PS, na.rm = TRUE),
    # total 
    total_avg = PN_avg + PS_avg, 
    # PN/PS
    ratio_avg = PN_avg/PS_avg, 
    .groups = "drop"
  ) %>%
  mutate(
    biofilm = factor(biofilm, levels = c("Methanogenic", "AGS")),
    exp_group = factor(exp_group, levels = c("Inner", "Outer", "Floccular", "S", "M", "L", "XL", "XXL"))
  )

summary_long <- summary_wide %>%
  pivot_longer(
    cols = c(PN_avg, PN_sd,
             PS_avg, PS_sd,
             total_avg, 
             ratio_avg),
    names_to = c("assay", ".value"),
    names_sep = "_"
  ) %>%
  mutate(
    y_label = if_else(assay == "ratio", "PN/PS", "\u00b5g/mgTSS"),
    y_label = factor(y_label, levels = c("\u00b5g/mgTSS", "PN/PS")),
    # Write out LB and TB
    extract = factor(extract, levels = c("TB", "LB")),
    extract = recode(extract, "LB" = "Loosely Bound", "TB" = "Tightly Bound"),
    # write out PN, PS, etc
    assay = factor(assay, levels = c("PS", "PN", "total", "ratio")),
    assay = recode(assay, "PN" = "Protein (PN)", "PS" = "Polysaccharide (PS)",
                   "total" = "Total EPS (PN + PS)", "ratio" = "PN/PS")
  ) 

# ------ Plot ------

max_tb_conc <- subset(summary_long, 
                      y_label == "\u00b5g/mgTSS" & extract == "Tightly Bound") %>%
  summarize(
    max_val = max(avg + sd, na.rm = TRUE)
    ) %>%
  pull(max_val)

max_lb_conc <- subset(summary_long, 
                      y_label == "\u00b5g/mgTSS" & extract == "Loosely Bound") %>%
  summarize(
    max_val = max(avg + sd, na.rm = TRUE)
  ) %>%
  pull(max_val)

max_tb_pnps <- subset(summary_long, 
                      y_label == "PN/PS" & extract == "Tightly Bound") %>%
  summarize(
    max_val = max(avg)
  ) %>%
  pull(max_val)

max_lb_pnps <- subset(summary_long, 
                      y_label == "PN/PS" & extract == "Loosely Bound") %>%
  summarize(
    max_val = max(avg)
  ) %>%
  pull(max_val)

  
p <- ggplot(summary_long, aes(x = exp_group, y = avg, fill = assay)) +
  
  # Concentration Plots
  geom_col(
    data = subset(summary_long, y_label == "\u00b5g/mgTSS"),
    position = "dodge",
    width = 0.8
  ) +
  geom_errorbar(
    data = subset(summary_long, y_label == "\u00b5g/mgTSS"),
    aes(ymin = avg - sd, ymax = avg + sd),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  
  # PN/PS plots
  geom_col(
    data = subset(summary_long, y_label == "PN/PS"),
    width = 0.5
  ) +
  geom_errorbar(
    data = subset(summary_long, y_label == "PN/PS"),
    aes(ymin = avg - sd, ymax = avg + sd),
    position = position_dodge(width = 0.5),
    width = 0.5/4
  ) +
  
  # Sizes
  # ggh4x::facet_grid2(
  ggh4x::facet_nested(
    extract + y_label ~ biofilm,
    scales = "free",
    independent = "y",
    switch = "y"
  ) +
  facetted_pos_scales(
    y = list(
      # TB
      scale_y_continuous(limits = c(0, 10)), # Methanogenic
      scale_y_continuous(limits = c(0, max_tb_conc)), # AGS
      scale_y_continuous(limits = c(0, max_tb_pnps), breaks = c(0, 2, 4)),
      scale_y_continuous(limits = c(0, max_tb_pnps), breaks = c(0, 2, 4)),
      # LB
      scale_y_continuous(limits = c(0, max_lb_conc), breaks = c(0, 3, 6, 9, 12)), # Methanogenic
      scale_y_continuous(limits = c(0, max_lb_conc), breaks = c(0, 3, 6, 9, 12)), # AGS
      scale_y_continuous(limits = c(0, max_lb_pnps), breaks = c(0, 2, 4)),
      scale_y_continuous(limits = c(0, max_lb_pnps), breaks = c(0, 2, 4))
    )
  ) +
  force_panelsizes(cols = c(0.35, 1), rows = c(1, 0.4, 1, 0.4)) +
  
  scale_fill_manual(
    values = c(
      "Polysaccharide (PS)" = "lightsalmon2",
      "Protein (PN)" = "lightblue",
      "Total EPS (PN + PS)" = "steelblue",
      "PN/PS" = "lightgray"
    )
  ) +
  
  labs(
    x = "Biofilm",
    y = NULL,
    fill = NULL
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1), 
    strip.placement = "outside",
    strip.background = element_blank()
  )


fname_out <- "./figures/Figure_2.tif"
ggsave(fname_out, plot = p, width = 6.5, height = 6, dpi = 300)
