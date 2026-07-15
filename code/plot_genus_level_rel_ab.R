rm(list = ls())
library(ComplexHeatmap)
library(MetBrewer)
source("./code/01_load_ps.R")

# number of rows to show
n_show <- 20

write2excel <- 0

# Cell height in inches (adjust as needed)
cell_h <- 0.2
cell_w <- 0.2 

# Font sizes
row_fontsize <- 10
col_fontsize <- 11

# Pseudo count to add when converting to log
# enables log(0)
pseudo <- 1e-6  

# ------ Define Data ------

# Relative Abundance
rel_df <- get_rel_wide(ps) %>%
  # Arrange taxa from largest to smallest abundance
  mutate(row_sum = rowSums(.)) %>%
  arrange(desc(row_sum)) %>%
  # Keep the top n_show
  head(., n = n_show) %>%
  dplyr::select(-row_sum) 

# sum of relative abundance in each sample
rel_sum <- colSums(rel_df)

# Convert to log
log_mat <- as.matrix(rel_df) %>%
  { log10(. + pseudo) }

# Metabolism
m_df <- rel_df %>%
  rownames_to_column(var = "Genus") %>%
  get_metabolism() 


#### ---- Plotting ------

# metabolism annotation
m_colors  <- c("P" = "#66C24A", "V" = "#EAEC3F")
m_annot <- rowAnnotation(
  df = m_df,
  # column names
  annotation_name_side = "bottom",
  annotation_name_rot = -60,
  # color
  col = col_list <- setNames(
    rep(list(m_colors), ncol(m_df)),
    colnames(m_df)
  ),
  na_col = NA, # no color for NA
  # legend
  show_legend = FALSE
)
# metabolism legend
lgd <- Legend(
  title = "Functional\nGroup",
  labels = c("Positive", "Variable"),
  legend_gp = gpar(fill = m_colors),
  nrow = 2,
  row_gap = unit(3, "mm")
)

## Relative Abundance

# Dimensions
n_cols <- ncol(log_mat)
n_rows <- nrow(log_mat)

# Labels
# Truncate names if they're longer than 25 characters (append with ...)
row_labels <- ifelse(
  nchar(rownames(log_mat)) > 25,
  paste0(substr(rownames(log_mat), 1, 25), "..."),  
  rownames(log_mat)
)
rownames(log_mat) = row_labels
# italicize classified genera
italic_rows <- !grepl("^(Unk|midas)", row_labels) 
# Apply
row_fontface <- ifelse(italic_rows, "italic", "plain")

# Legend Colors
ht_colors <- met.brewer(taxa_pal, type = "continuous")
# Display legend ticks
break_values <- c(0, 0.1, 1, 10, 73) # %
breaks_log_display <- log10(break_values + pseudo) # Log (%)
# Add % symbol to top break
breaks_rel_display <- replace(
  as.character(break_values),
  length(break_values),
  paste0(tail(break_values, 1), "%")
)

ht <- Heatmap(
  log_mat,
  # columns
  show_column_names = TRUE,  
  column_names_rot = -60,
  column_title = NULL, #"Relative Abundance",
  cluster_columns = FALSE, # changes sample order
  # heatmap legend
  col = ht_colors, 
  show_heatmap_legend = TRUE, 
  heatmap_legend_param = list(
    at = breaks_log_display,
    labels = breaks_rel_display,
    title = "Relative Abundance", 
    title_position = "leftcenter-rot",
    legend_height = unit(7, "cm")
  ),
  # Annotations
  right_annotation = m_annot,
  # Display size
  width  = unit(n_cols * cell_w, "inches"),
  height = unit(n_rows * cell_h, "inches"),
  row_names_gp = gpar(fontsize = row_fontsize, fontface = row_fontface),
  column_names_gp = gpar(fontsize = col_fontsize)
)

# Figure output location
fname_rel <- "./figures/genus_level_rel_ab.png"

# Draw combined heatmap
png(fname_rel,
    width = 4,  # width in inches; can adjust
    height = 5, # height in inches; can adjust
    units = "in", res = 300)
draw(ht, heatmap_legend_side = "left") 
# metabolism legend
draw(lgd, x = unit(0.03, "npc"), y = unit(0.25, "npc"), just = c("left", "top")) 
dev.off()

## Check what percent of relative abundance is included in plot
message(paste("Heatmap Min:", round(min(rel_sum), 2), "%"))
message(paste("Heatmap Max:", round(max(rel_sum), 2), "%"))
