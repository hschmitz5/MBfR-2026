rm(list = ls())
library(tidyverse)
library(readxl)
library(writexl)
source("./code/01_load_ps.R")

# define sample names
group_levels <- c("Inner", "Outer", "Floccular", "S", "M", "L", "XL", "XXL")

raw_df <- read_excel("./data/Rheometry_MBfR_AGS.xlsx", sheet = "input", skip = 1)

modulus <- raw_df %>%
  select(-freq_hz) %>%
  pivot_longer(
    cols = c(G_1:G_7, G2_1:G2_7),
    names_to = c(".value", "replicate"),
    names_pattern = "(G2?)_(\\d)"
  ) %>%
  filter(
    freq_rad < 40,
    !is.na(G)
    ) %>%
  # * indicates apparent groups (with voids)
  mutate(
    exp_group = factor(exp_group, levels = group_levels),
    exp_group = recode(exp_group, 
                       "S" = "S*", "M" = "M*", "L" = "L*", "XL" = "XL*", "XXL" = "XXL*")
  )

# ------ ANOVA ------

modulus_sub <- modulus %>%
  filter(exp_group %in% c("Inner", "Outer", "Floccular")) 

# storage modulus
G_res <- aov(modulus_sub$G ~ modulus_sub$exp_group)

G_res_pairwise <- TukeyHSD(G_res)$`modulus_sub$exp_group` %>%
  as.data.frame() %>%
  rownames_to_column(var = "Comparison") %>%
  mutate(mod_type = "G", .before = 1)

# loss modulus
G2_res <- aov(modulus_sub$G2 ~ modulus_sub$exp_group)

G2_res_pairwise <- TukeyHSD(G2_res)$`modulus_sub$exp_group` %>%
  as.data.frame() %>%
  rownames_to_column(var = "Comparison") %>%
  mutate(mod_type = "G2", .before = 1)

res_pair <- bind_rows(G_res_pairwise, G2_res_pairwise)

fname_out <- "./data/G_stats.xlsx"
write_xlsx(res_pair, path = fname_out)

# ------ Summarize mean across replicates for plotting ------

mod_summary_wide <- modulus %>%
  group_by(biofilm, exp_group, freq_rad) %>%
  summarize(
    G_avg = mean(G),
    G_sd = sd(G),
    G2_avg = mean(G2),
    G2_sd = sd(G2),
    .groups = "drop"
  ) 

mod_summary_long <- mod_summary_wide %>%
  pivot_longer(
    cols = c(G_avg, G_sd, G2_avg, G2_sd),
    names_to = c("measure", ".value"),
    names_pattern = "(G2?|G2?)_(avg|sd)"
  ) %>%
  mutate(
    # convert units to kPa (originally in Pa)
    avg = avg/1000, 
    sd = sd/1000,
    # change display names and order
    measure = factor(measure, levels = c("G", "G2")),
    measure = recode(measure,"G"="Storage Modulus (G')","G2"='Loss Modulus (G")')
  )


#### Plot

p1 <- ggplot(mod_summary_long, aes(x = freq_rad, y = avg, color = exp_group)) +
  geom_point() +
  geom_line(aes(group = exp_group)) +
  geom_errorbar(
    aes(ymin = pmax(avg - sd, 0), ymax = avg + sd),
    width = 0.2
  ) +
  facet_wrap(~measure, scales = "free_y", nrow = 1) +
  scale_y_log10() +  # changes scale distribution, not values
  scale_color_manual(
    name = "Fraction",
    values = c("darkseagreen3", "steelblue", "gray", met.brewer("Greek", 5))
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

fname_out <- "./figures/Figure_1.tif"
ggsave(fname_out, plot = p1, width = 6.5, height = 2.5, dpi = 300)
