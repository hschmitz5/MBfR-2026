rm(list = ls())
library(tidyverse)
library(patchwork)

# File names for concentration data
fname_pn_mbfr    <- paste0("./data/EPS/PN_conc_mbfr.rds") 
fname_polys_mbfr <- paste0("./data/EPS/PS_conc_mbfr.rds") 
fname_pn_ags     <- paste0("./data/EPS/PN_conc_ags.rds") 
fname_polys_ags  <- paste0("./data/EPS/PS_conc_ags.rds") 
  
# Calculate average and std of replicates
group_data <- function(fname) {
  df <- readRDS(fname) %>%
    rename(.,"biofilm" = if (hasName(., "region")) "region" else "size") %>%
    group_by(biofilm, extract) %>% # region
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
  PN %>% select(biofilm, extract, PN_avg = avg, PN_sd = sd), 
  PS %>% select(biofilm, extract, PS_avg = avg, PS_sd = sd), 
  by = c("biofilm", "extract")
  ) %>%
  mutate(
    total = PN_avg + PS_avg,
    PNPS = PN_avg/PS_avg,
    sd = NA
    ) 

# Combine into single data frame
df <- bind_rows(
  'Protein (PN)' = PN,
  'Polysaccharide (PS)' = PS,
  'Total EPS (PN + PS)' = df_wide %>% select(biofilm, extract, avg = total, sd),
  .id = "assay"
  ) %>%
  mutate(
    assay = factor(assay, levels = c("Polysaccharide (PS)", "Protein (PN)", "Total EPS (PN + PS)"))
    ) 

# Calculate PN/PS
PNPS <- df_wide %>% 
  select(biofilm, extract, avg = PNPS) 


# ------ Plot ------

# Determine maximum avg + sd
max_y1 <- df %>%
  filter(assay != "Total EPS (PN + PS)") %>%
  summarise(max_y = max(avg + sd)) %>%
  pull(max_y)

max_y2 <- df %>%
  filter(assay == "Total EPS (PN + PS)") %>%
  summarise(max_y = max(avg)) %>%
  pull(max_y)

max_y <- ceiling(
  max(max_y1, max_y2)
  )

p <- ggplot(data = df, aes(x = biofilm, y = avg, fill = assay)) +
  geom_col(position = "dodge", width = 0.8) +
  geom_errorbar(
    aes(ymin = avg - sd, ymax = avg + sd),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  facet_wrap(~extract, scales = "free_y") + 
  # ylim(0, max_y) +
  labs(
    y = expression(paste(mu, "g/mgTSS")),
    x = NULL,
    fill = NULL # legend titles
  ) +
  scale_fill_manual(
    values = c(
      "Polysaccharide (PS)" = "lightsalmon2",
      "Protein (PN)"        = "lightblue",
      "Total EPS (PN + PS)"      = "steelblue" 
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.background = element_rect(
      colour = NA # facet label outline
    ),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) 

annot <- ggplot(data = PNPS, aes(x = biofilm, y = avg, fill = "PN/PS")) +
  geom_col(position = "dodge", width = 0.5) +
  facet_wrap(~extract) +
  scale_y_continuous(
    breaks = c(0, 2, 4)
    ) +
  labs(
    x = "Biofilm",
    y = NULL, 
    fill = NULL
    ) +
  scale_fill_manual(
    values = "lightgray",
    labels = expression(frac(PN, PS))
  ) +
  theme_classic(base_size = 12) +
  theme(
    strip.text  = element_blank(),
    legend.justification = "left",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p2 <- p / annot +
  plot_layout(heights = c(4, 1.5))

fname_out <- "./figures/EPS_comparison.png"
ggsave(fname_out, plot = p2, width = 6.5, height = 3, dpi = 300)
