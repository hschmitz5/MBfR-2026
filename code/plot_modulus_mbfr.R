rm(list = ls())
library(readxl)
library(tidyverse)
library(patchwork)
library(MetBrewer)

fname_in  <- "./data/RheometryApr142026.xlsx"
fname_out <- "./figures/moduli.png"

modulus <- read_excel(fname_in, sheet = "input", skip = 1) %>%
  mutate(
    region = recode(region,"B" = "inner","M" = "outer")
  )

#### Plot

p1 <- ggplot(modulus, aes(x = freq_rad, y = G_avg, color = region)) +
  geom_point() +
  geom_line(aes(group = region)) +
  geom_errorbar(
    aes(ymin = pmax(G_avg - G_sd, 0), ymax = G_avg + G_sd),
    width = 0.2
  ) +
  scale_color_manual(
    values = c("inner" = "darkseagreen", "outer" = "lightsalmon3")
  ) +
  labs(
    x = "Frequency [rad/s]",
    y = "Storage Modulus [Pa]",
    color = "Region"
  ) 

p2 <- ggplot(modulus, aes(x = freq_rad, y = G2_avg, color = region)) +
  geom_point() +
  geom_line(aes(group = region)) +
  geom_errorbar(
    aes(ymin = pmax(G2_avg - G2_sd, 0), ymax = G2_avg + G2_sd),
    width = 0.2
  ) +
  scale_color_manual(
    values = c("inner" = "darkseagreen", "outer" = "lightsalmon3")
  ) +
  labs(
    x = "Frequency [rad/s]",
    y = "Loss Modulus [Pa]",
    color = "Region"
  ) 


# horizontal
p <- p1 + p2 + 
  plot_layout(guides = "collect") & 
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom") &
  guides(
    color = guide_legend(
      title.position = "bottom",
      title.hjust = 0.5 # centers title
      )
  )

ggsave(fname_out, plot = p, width = 8, height = 3, dpi = 600)
