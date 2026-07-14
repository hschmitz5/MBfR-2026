rm(list = ls())
library(tidyverse)
library(patchwork)
source("./code/01_load_ps.R")

# define sample names
region <- c("inner", "outer")

fname_in <- "./data/RheometryApr142026.xlsx"
modulus <- read_excel(fname_in, sheet = "input", skip = 1) %>%
  select(region, freq_rad, G_avg, G_sd, G2_avg, G2_sd) %>%
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
    measure = factor(measure, levels = c("G", "G2")),
    measure = recode(measure,"G"="Storage Modulus","G2"="Loss Modulus")
    )

modulus_subset <- modulus %>%
  filter(freq_rad == 0.1) %>%
  select(-freq_rad) 

#### Plot

p <- ggplot(modulus, aes(x = freq_rad, y = avg, color = region)) +
  geom_point() +
  geom_line(aes(group = region)) +
  geom_errorbar(
    aes(ymin = pmax(avg - sd, 0), ymax = avg + sd),
    width = 0.2
  ) +
  facet_wrap(~measure, scales = "free_y", nrow = 1) +
  scale_color_manual(
    name = "Region",
    values = c("inner" = "darkseagreen", "outer" = "lightsalmon3")
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

fname_out <- "./figures/moduli.png"
ggsave(fname_out, plot = p, width = 6.5, height = 2.25, dpi = 300)

p_sub <- ggplot(modulus_subset, aes(x = region, y = avg, fill = measure)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_errorbar(
    aes(ymin = avg - sd, ymax = avg + sd),
    width = 0.2,
    position = position_dodge(width = 0.6)
  ) +
  labs(
    title = "Frequency = 0.1 rad/s",
    x = "Region",
    y = "Modulus (kPa)"
  ) +
  scale_fill_manual(
    values = c("plum4", "lightgray")
  ) +
  theme_classic(base_size = 12) +
  theme(
    legend.title = element_blank()
  )

fname_out <- "./figures/moduli_subset.png"
ggsave(fname_out, plot = p_sub, width = 5, height = 2.25, dpi = 300)
