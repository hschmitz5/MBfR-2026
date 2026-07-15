rm(list = ls())
library(tidyverse)
library(ggh4x)

# File names for concentration data
fname_pn_mbfr    <- paste0("./data/EPS/PN_conc_mbfr.rds") 
fname_polys_mbfr <- paste0("./data/EPS/PS_conc_mbfr.rds") 

# Calculate average and std of replicates
group_data <- function(fname) {
  df <- readRDS(fname) %>%
    group_by(extract, region) %>% # region
    summarize(
      avg = mean(C_TSS),
      sd = sd(C_TSS),
      .groups = "drop"
    ) %>%
    mutate(
      region = factor(region, levels = c("inner", "outer")),
      extract = recode(extract,"LB" = "Loosely Bound","TB" = "Tightly Bound"),
      extract = factor(extract, levels = c("Tightly Bound", "Loosely Bound"))
    )
}
# Apply function to each assay
PN <- group_data(fname_pn_mbfr)
PS <- group_data(fname_polys_mbfr)

# Calculate PN + PS and PN/PS
df_wide <- left_join(
  PN %>% select(extract, region, PN_avg = avg, PN_sd = sd), 
  PS %>% select(extract, region, PS_avg = avg, PS_sd = sd), 
  by = c("extract", "region")
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
  'Total EPS (PN + PS)' = df_wide %>% select(extract, region, avg = total, sd),
  .id = "assay"
  ) %>%
  mutate(plot_type = "\u00b5g/mgTSS") %>%
  select(plot_type, assay, extract, region, avg, sd) 

# Calculate PN/PS
PNPS <- df_wide %>% 
  mutate(
    plot_type = "PN/PS",
    assay = "PN/PS",
    sd = NA
  ) %>%
  select(plot_type, assay, extract, region, avg = PNPS, sd) 
  
df_all <- bind_rows(df_conc, PNPS) %>%
  mutate(
    assay = factor(assay, levels = c("Polysaccharide (PS)", "Protein (PN)", "Total EPS (PN + PS)", "PN/PS"))
    )

# ------ Plot ------


p <- ggplot(df_all, aes(x = region, y = avg, fill = assay)) +
  
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
  facet_grid(
    plot_type ~ extract,
    scales = "free_y",
    switch = "y"
  ) +
  facetted_pos_scales(
    y = list(
      scale_y_continuous(),                   # Concentration row
      scale_y_continuous(breaks = c(0, 2.5, 5)) # PN/PS row
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
    x = "Region",
    y = NULL,
    fill = NULL
  ) +
  
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.placement = "outside",
    strip.background = element_blank()
  )


fname_out <- "./figures/EPS.png"
ggsave(fname_out, plot = p, width = 6.5, height = 3, dpi = 300)
