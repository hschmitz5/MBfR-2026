rm(list = ls())
library(tidyverse)
library(readxl)
library(cowplot)
source("./code/01_load_ps.R")

# define sample names
biofilm_levels <- c("inner", "outer", "XS")

fname_in <- "./data/RheometryApr142026.xlsx"
modulus <- read_excel(fname_in, sheet = "input", skip = 1) %>%
  select(biofilm, freq_rad, G_avg, G_sd, G2_avg, G2_sd) %>%
  pivot_longer(
    cols = c(G_avg, G_sd, G2_avg, G2_sd),
    names_to = c("measure", ".value"),
    names_pattern = "(G2?|G2?)_(avg|sd)"
  ) %>%
  mutate(
    # convert units to MPa (originally in Pa)
    avg = avg/1e3, 
    sd = sd/1e3,
    # change display names and order
    biofilm = factor(biofilm, levels = biofilm_levels),
    measure = factor(measure, levels = c("G", "G2")),
    measure = recode(measure,"G"="Storage Modulus","G2"="Loss Modulus")
    ) %>%
  filter(freq_rad < 110)

modulus_subset <- modulus %>%
  filter(freq_rad == 0.1) %>%
  select(-freq_rad) 

#### Plot

p1 <- ggplot(modulus, aes(x = freq_rad, y = avg, color = biofilm)) +
  geom_point() +
  geom_line(aes(group = biofilm)) +
  geom_errorbar(
    aes(ymin = pmax(avg - sd, 0), ymax = avg + sd),
    width = 0.2
  ) +
  facet_wrap(~measure, scales = "free_y", nrow = 1) +
  scale_color_manual(
    name = "Biofilm",
    values = c("inner" = "steelblue", "outer" = "lightsalmon2", "XS" = "black")
  ) +
  labs(
    x = "Frequency (rad/s)",
    y = "Modulus (kPa)",
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "right",
    strip.background = element_rect(
      colour = NA # facet label outline
      )
    )

p2 <- ggplot(modulus_subset, aes(x = biofilm, y = avg, fill = measure)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_errorbar(
    aes(ymin = avg - sd, ymax = avg + sd),
    width = 0.2,
    position = position_dodge(width = 0.6)
  ) +
  labs(
    title = "Frequency = 0.1 rad/s",
    x = "Biofilm",
    y = "Modulus (kPa)"
  ) +
  scale_fill_manual(
    values = c("plum4", "lightgray")
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.title = element_blank()
  )

# arrange two plots into one column
p <- plot_grid(
  p1, p2,
  labels = "AUTO", ncol = 1, rel_widths = c(6.5, 5)
)

fname_out <- "./figures/moduli.png"
ggsave(fname_out, plot = p, width = 6.5, height = 5, dpi = 300)

# fname_out <- "./figures/moduli.png"
# ggsave(fname_out, plot = p1, width = 6.5, height = 2.25, dpi = 300)

# fname_out <- "./figures/moduli_subset.png"
# ggsave(fname_out, plot = p2, width = 5, height = 2.25, dpi = 300)
