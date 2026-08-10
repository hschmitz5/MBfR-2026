rm(list = ls())
library(tidyverse)

import_data <- function(PN_fname, PS_fname, biofilm_type, group_name) {
  PN <- readRDS(PN_fname) %>%
    select(extract, group_name, replicate, PN = C_TSS) 
  
  PS <- readRDS(PS_fname) %>%
    select(extract, group_name, replicate, PS = C_TSS)
  
  df_joined <- left_join(PN, PS, by = c("extract", group_name, "replicate")) %>%
    mutate(biofilm = biofilm_type) %>%
    select(biofilm, extract, exp_group = group_name, replicate, PN, PS)
}

ags_data <- import_data("./data/EPS/PN_conc_ags.rds", "./data/EPS/PS_conc_ags.rds", "AGS", "size") 

mbfr_data <- import_data("./data/EPS/PN_conc_mbfr.rds", "./data/EPS/PS_conc_mbfr.rds", "MBfR", "region") %>%
  mutate(replicate = as.character(replicate))

# Combine both data sets
eps <- bind_rows(ags_data, mbfr_data)


# ------ t test ------

eps_var = c("PN", "PS") 

res_overall <- map_dfr(c("LB", "TB"), \(ext) {
  df <- filter(eps, extract == ext)
  
  tibble(
    extract = ext,
    variable = eps_var,
    p.value = map_dbl(
      eps_var,
      \(var) t.test(df[[var]] ~ df$biofilm)$p.value
    )
  ) %>%
    mutate(p.adjusted = p.adjust(p.value, method = "BH"))
})


# ------ pairwise ANOVA ------

exp_group_levels <- c("Inner", "Outer", "Floccular", "S", "M", "L", "XL", "XXL")

res_pairwise_long <- map_dfr(c("LB", "TB"), \(ext) {
  
  df <- filter(eps, extract == ext) %>%
    mutate(
      exp_group = factor(exp_group, levels = exp_group_levels)
      )
  
  map_dfr(eps_var, \(var) {
    
    aov_res <- aov(df[[var]] ~ df$exp_group)
    
    TukeyHSD(aov_res)[[1]] %>%
      as.data.frame() %>%
      tibble::rownames_to_column("comparison") %>%
      tidyr::separate(
        comparison,
        into = c("group_1", "group_2"),
        sep = "-"
        ) %>%
      rename(p.adj = `p adj`) %>%
      mutate(
        extract = ext,
        variable = var,
        .before = 1
        )
    })
  }) %>%
  mutate(
    group_1 = factor(group_1, levels = exp_group_levels),
    group_2 = factor(group_2, levels = exp_group_levels)
    )

res_pairwise <- map(c("LB", "TB"), \(ext) {
  map(eps_var, \(var) {
    res_pairwise_long %>%
      filter(
        extract == ext,
        variable == var
        ) %>%
      select(group_1, group_2, p.adj) %>%
      pivot_wider(
        names_from = group_1,
        values_from = p.adj
        )
    }) %>%
    set_names(eps_var)
  }) %>%
  set_names(c("LB", "TB"))

res_pairwise$LB$PN
res_pairwise$LB$PS

res_pairwise$TB$PN
res_pairwise$TB$PS
