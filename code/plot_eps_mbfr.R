rm(list = ls())
library(tidyverse)
library(patchwork)

root_name   <- "fresh"

# -----------------------

# Format Data

fname_pn    <- paste0("./data/",root_name,"_PN_conc.rds")
fname_polys <- paste0("./data/",root_name,"_PS_conc.rds")

# Calculate average and std of replicates
group_data <- function(fname) {
  df <- readRDS(fname) %>%
    group_by(biofilm, extract) %>%
    summarize(
      avg = mean(C_TSS),
      sd = sd(C_TSS),
      .groups = "drop"
    )
}
# Apply function to each assay
PN <- group_data(fname_pn) 
PS <- group_data(fname_polys)

# Combine into single data frame
df <- bind_rows(
    PN = PN,
    PS = PS,
    .id = "assay"
  ) %>%
  mutate(
    extract = recode(extract,"LB" = "Loosely Bound","TB" = "Tightly Bound")
  )

# Calculate PN/PS and total EPS
df_wide <- df %>%
  pivot_wider(
    names_from = assay,
    values_from = c(avg, sd),
    names_glue = "{assay}_{.value}"
  ) %>%
  mutate(
    PNPS = PN_avg/PS_avg,
    total_avg = PN_avg + PS_avg
  )

# Define total EPS
tot <- df_wide %>%
  transmute(assay = "total", biofilm, extract, avg = total_avg)

# Determine maximum avg + sd
max_y <- ceiling(
  max(df$avg + df$sd)
)

# Get rid of acronyms
df <- df %>%
  mutate(
    assay = recode(assay, "PN" = "protein (PN)", "PS" = "polysaccharide (PS)")
  )

# -------------------------------

# Make Plot
p <- ggplot(data = df, aes(x = biofilm, y = avg, fill = assay)) +
  geom_col(position = "dodge", width = 0.6) +
  geom_errorbar(
    aes(ymin = avg - sd, ymax = avg + sd),
    position = position_dodge(width = 0.6),
    width = 0.2
  ) +
  geom_point(data = tot, aes(x = biofilm, y = avg, color = assay), shape = 20) +
  facet_wrap(~extract, nrow=1, strip.position = "top") +
  ylim(0, max_y) +
  labs(
    x = "Region",
    y = expression(paste(mu, "g/mgTSS")),
    color = NULL,
    fill = NULL
  ) +
  scale_fill_manual(
    values = c(
      "protein (PN)" = "gray",
      "polysaccharide (PS)" = "lightblue"
    )
  ) +
  scale_color_manual(
    values = c("total" = "black")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    strip.background = element_rect(
      fill = "white",
      colour = "lightgray"
    )
  ) +
  guides(
    fill = guide_legend(order = 1),
    color = guide_legend(order = 2)
  ) 

annot <- ggplot(data = df_wide) +
  geom_text(aes(x = biofilm, y = "PN/PS", label = round(PNPS,1))) +
  scale_fill_gradient(low="white", high="lightgray") +
  facet_wrap(~extract, nrow=1, strip.position = "top") +
  theme(
    legend.position = "none",
    panel.grid  = element_blank(),
    axis.title  = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks  = element_blank(),
    strip.text  = element_blank()
  )

p2 <- p / annot +
  plot_layout(heights = c(4, 1)) 

fname_out <- paste0("./figures/",root_name,"_EPS.png")
ggsave(fname_out, plot = p2, width = 8, height = 4, dpi = 600)
