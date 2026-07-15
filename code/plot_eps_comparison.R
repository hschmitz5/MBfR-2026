rm(list = ls())
library(tidyverse)
library(ggh4x)

# File names for concentration data
fname_pn_mbfr    <- paste0("./data/EPS/PN_conc_mbfr.rds") 
fname_polys_mbfr <- paste0("./data/EPS/PS_conc_mbfr.rds") 
fname_pn_ags     <- paste0("./data/EPS/PN_conc_ags.rds") 
fname_polys_ags  <- paste0("./data/EPS/PS_conc_ags.rds") 
  
# Calculate average and std of replicates
group_data <- function(fname) {
  df <- readRDS(fname) %>%
    rename(.,"biofilm" = if (hasName(., "region")) "region" else "size") %>%
    group_by(extract, biofilm) %>% # region
    summarize(
      avg = mean(C_TSS),
      sd = sd(C_TSS),
      .groups = "drop"
    ) %>%
    mutate(
      biofilm = factor(biofilm, levels = c("inner", "outer", "XS", "S", "M", "L", "XL", "XXL")),
      extract = recode(extract,"LB" = "Loosely Bound","TB" = "Tightly Bound"),
      extract = factor(extract, levels = c("Tightly Bound", "Loosely Bound"))
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
  PN %>% select(extract, biofilm, PN_avg = avg, PN_sd = sd), 
  PS %>% select(extract, biofilm, PS_avg = avg, PS_sd = sd), 
  by = c("extract", "biofilm")
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
  'Total EPS (PN + PS)' = df_wide %>% select(extract, biofilm, avg = total, sd),
  .id = "assay"
) %>%
  mutate(plot_type = "\u00b5g/mgTSS") %>%
  select(plot_type, assay, extract, biofilm, avg, sd) 

# Calculate PN/PS
PNPS <- df_wide %>% 
  mutate(
    plot_type = "PN/PS",
    assay = "PN/PS",
    sd = NA
  ) %>%
  select(plot_type, assay, extract, biofilm, avg = PNPS, sd) 

df_all <- bind_rows(df_conc, PNPS) %>%
  mutate(
    assay = factor(assay, levels = c("Polysaccharide (PS)", "Protein (PN)", "Total EPS (PN + PS)", "PN/PS"))
  )

# ------ Plot ------


p <- ggplot(df_all, aes(x = biofilm, y = avg, fill = assay)) +
  
  # Concentration Plots
  geom_col(
    data = subset(df_all, plot_type == "\u00b5g/mgTSS"),
    position = "dodge",
    width = 0.8
  ) +
  geom_errorbar(
    data = subset(df_all, plot_type == "\u00b5g/mgTSS"),
    aes(ymin = avg - sd, ymax = avg + sd),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  
  # PN/PS plots
  geom_col(
    data = subset(df_all, plot_type == "PN/PS"),
    width = 0.5
  ) +
  
  # Sizes
  ggh4x::facet_grid2(
    plot_type ~ extract,
    scales = "free",
    switch = "y",
    independent = "y"
  ) +
  facetted_pos_scales(
    y = list(
      scale_y_continuous(),   
      scale_y_continuous(), 
      scale_y_continuous(breaks = c(0, 2.5, 5)), # PN/PS row
      scale_y_continuous(breaks = c(0, 2.5, 5))  # PN/PS row
    )
  ) +
  force_panelsizes(rows = c(1, 1/3), cols = c(1, 1)) +
  
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
ggsave(fname_out, plot = p, width = 6.5, height = 3, dpi = 300)
