# NOTE: model forces intercept to intersect (0,0)

# Verify quality of predicted output 
# by verifying that the data is within the working range,
# and check the R-squared in the model summary

rm(list = ls())
library(tidyverse)
library(readxl)
library(ggplot2)

# change to process each sheet
sheet_name <- "PN"  # protein (PN) or polysaccharide (PS)

# File name for absorbance data
fname_in   <- "./data/EPS/EPS_ctr_fresh.xlsx"

# Used to convert concentration to include TSS/VSS
extract_volume = 10 # mL

# ----------------------------------------


# Read in absorbance data
df <- read_excel(fname_in, sheet = sheet_name, skip = 1) %>%
  mutate(
    region = recode(region,"B" = "Inner","M" = "Outer")
  ) 

# define standards
std <- df %>%
  filter(is.finite(C)) %>%
  select(sample,A,C)

# define samples
sam <- df %>%
  filter(is.na(C)) %>%
  select(-C)

# set poly_degree based on sheet_name
poly_degree <- switch(sheet_name,
                      "PN" = 3,
                      "PS" = 2,
                      stop("invalid sheet_name"))

## fit concentration as a function of absorbance
# 0 forces the equation to intercept (0,0)
model <- lm(C ~ 0 + poly(A, poly_degree, raw = TRUE), data = std)
msum <- summary(model)

# predict concentrations
sam$C0 <- predict(model, newdata = sam)
sam$C_TSS <- sam$C0*extract_volume/sam$TSS
#sam$C_VSS <- sam$C0*extract_volume/sam$VSS

# save the sample data
saveRDS(sam, file = paste0("./data/EPS/",sheet_name,"_conc_mbfr.rds"))


#### Plot Fit Data

# fit line for plotting
A_seq <- seq(min(std$A), max(std$A), length.out = 200)
fit_df <- data.frame(
  A = A_seq,
  C = predict(model, newdata = data.frame(A = A_seq))
)

# plot
ggplot() +
  geom_line(data = fit_df, aes(C, A, color = "Fit"), linewidth = 1) +
  geom_point(data = std, aes(C, A, color = "Standards")) +
  geom_point(data = sam, aes(C0, A, color = "Samples"), shape = 3) +
  scale_color_manual(
    values = c(
      "Fit" = "lightblue",
      "Standards" = "black",
      "Samples" = "red"
    ),
    labels = c(
      "Fit" = paste0("Fit (R² = ", round(msum$r.squared, 4), ")"),
      "Standards" = "Standards",
      "Samples" = "Samples"
    )
  ) +
  labs(
    color = NULL,
    x = bquote("Concentration [" * mu * "g" * .(sheet_name) * "/mL]"),
    y = "Absorbance"
  ) +
  theme_minimal(base_size = 12) +
  theme(aspect.ratio = 0.7)

fit_plot <- paste0("./figures/EPS/",sheet_name,"_fit.png")
ggsave(fit_plot, height = 2.5, width = 6, dpi = 600)