rm(list = ls())
library(tidyverse)
library(ggh4x)

# File names for concentration data
fname_pn_mbfr    <- paste0("./data/EPS/PN_conc_mbfr.rds") 
fname_polys_mbfr <- paste0("./data/EPS/PS_conc_mbfr.rds") 
fname_pn_ags     <- paste0("./data/EPS/PN_conc_ags.rds") 
fname_polys_ags  <- paste0("./data/EPS/PS_conc_ags.rds") 

biofilm <- data.frame(
  exp_category = c("Inner", "Outer", "Floccular", "S", "M", "L", "XL", "XXL"),
  type = c("MBfR", "MBfR", "AGS", "AGS", "AGS", "AGS", "AGS", "AGS")
)

# Calculate average and std of replicates
group_data <- function(fname) {
  df <- readRDS(fname) %>%
    rename(.,"exp_category" = if (hasName(., "region")) "region" else "size") %>%
    group_by(extract, exp_category) %>% # region
    summarize(
      avg = mean(C_TSS),
      sd = sd(C_TSS),
      .groups = "drop"
    ) %>%
    mutate(
      exp_category = factor(exp_category, levels = c("Inner", "Outer", "Floccular", "S", "M", "L", "XL", "XXL")),
      extract = recode(extract,"LB" = "Loosely Bound","TB" = "Tightly Bound"),
      extract = factor(extract, levels = c("Tightly Bound", "Loosely Bound")),
      biofilm = biofilm$type[as.numeric(exp_category)],
      biofilm = factor(biofilm, levels = c("MBfR", "AGS"))
    )
}
# Apply function to each assay
PN <- bind_rows(
  group_data(fname_pn_mbfr), 
  group_data(fname_pn_ags) 
)
PS <- bind_rows(
  group_data(fname_polys_mbfr),
  group_data(fname_polys_ags)
)

# Calculate PN + PS and PN/PS
df_wide <- left_join(
  PN %>% select(extract, biofilm, exp_category, PN_avg = avg, PN_sd = sd), 
  PS %>% select(extract, biofilm, exp_category, PS_avg = avg, PS_sd = sd), 
  by = c("extract", "biofilm", "exp_category")
) %>%
  mutate(
    total = PN_avg + PS_avg,
    PNPS = PN_avg/PS_avg,
    sd = NA
  ) 

# Combine into single data frame
df_conc <- bind_rows(
  'Protein (PN)' = PN,
  'Polysaccharide (PS)' = PS,
  'Total EPS (PN + PS)' = df_wide %>% select(extract, biofilm, exp_category, avg = total, sd),
  .id = "assay"
) %>%
  mutate(y_label = "\u00b5g/mgTSS") %>%
  select(y_label, assay, extract, biofilm, exp_category, avg, sd) 

# Calculate PN/PS
PNPS <- df_wide %>% 
  mutate(
    y_label = "",
    assay = "PN/PS",
    sd = NA
  ) %>%
  select(y_label, assay, extract, biofilm, exp_category, avg = PNPS, sd) 

df_all <- bind_rows(df_conc, PNPS) %>%
  mutate(
    y_label = factor(y_label, levels = c("\u00b5g/mgTSS", "")),
    assay = factor(assay, levels = c("Polysaccharide (PS)", "Protein (PN)", "Total EPS (PN + PS)", "PN/PS"))
  )

# ------ Plot ------

p <- ggplot(df_all, aes(x = exp_category, y = avg, fill = assay)) +
  
  # Concentration Plots
  geom_col(
    data = subset(df_all, y_label == "\u00b5g/mgTSS"),
    position = "dodge",
    width = 0.8
  ) +
  geom_errorbar(
    data = subset(df_all, y_label == "\u00b5g/mgTSS"),
    aes(ymin = avg - sd, ymax = avg + sd),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  
  # PN/PS plots
  geom_col(
    data = subset(df_all, y_label == ""),
    width = 0.5
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
      scale_y_continuous(limits = c(0, 80)), # TB MBfR
      scale_y_continuous(limits = c(0, 80)), # TB AGS
      scale_y_continuous(breaks = c(0, 2.5, 5)),
      scale_y_continuous(limits = c(0, 5), breaks = c(0, 2.5, 5)),
      # LB
      scale_y_continuous(limits = c(0, 15)), # LB MBfR
      scale_y_continuous(limits = c(0, 15)), # LB AGS
      scale_y_continuous(breaks = c(0, 2.5, 5)),
      scale_y_continuous(limits = c(0, 5), breaks = c(0, 2.5, 5))
    )
  ) +
  force_panelsizes(cols = c(1/3, 1), rows = c(1, 0.4, 1, 0.4)) +
  
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


fname_out <- "./figures/EPS_comparison.png"
ggsave(fname_out, plot = p, width = 6.5, height = 6, dpi = 300)
