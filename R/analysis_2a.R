# analysis 

library(dplyr)
library(tidyverse)
library(vegan)
library(mobr)
library(ggplot2)
library(betapart)
library(lme4)
library(performance)
library(tibble)
library(scales)
library(DHARMa)
library(see)
library(glmmTMB)

# DATA
wet_all_summer <- readRDS("Intermediate_data/all_wet_check_dominants_removed_summer.RDS")

wet_all_fall <- readRDS("Intermediate_data/all_wet_check_dominants_removed_fall.RDS")

wet_all_spring <- readRDS("Intermediate_data/all_wet_check_dominants_removed_spring.RDS")

wet_all_winter <- readRDS("Intermediate_data/all_wet_check_dominants_removed_winter.RDS")

env_summer <- readRDS("Data/earth_engine_env_data/env_matrix.RDS")

# need spatial coords for the LOCALITY_IDs
wet_coords <- readRDS("Intermediate_data/locality_ids_long_lat.RDS")

env <- readRDS("Data/earth_engine_env_data/env_matrix.RDS")

env <- env %>%
  left_join(wet_coords, by = "LOCALITY_ID")

## Creating a dataframe with just the LOCALITY_ID and ECOREGION
## to append to the wet_all df

ecoregion_loc <- env %>%
  select(LOCALITY_ID, NA_L1NAME)

wet_all_summer <- wet_all_summer %>%
  inner_join(ecoregion_loc, by = "LOCALITY_ID")

wet_all_fall <- wet_all_fall %>%
  inner_join(ecoregion_loc, by = "LOCALITY_ID")

wet_all_spring <- wet_all_spring %>%
  inner_join(ecoregion_loc, by = "LOCALITY_ID")

wet_all_winter <- wet_all_winter %>%
  inner_join(ecoregion_loc, by = "LOCALITY_ID") 

######### OBJECTIVE 3 - GIW scale ############
######### Within GIW beta diversity analysis ##########

# We just need to slightly alter the function, to calculate betas WITHIN a single wetland
# to see if we have high turnover simply between checklists at the same wetland


calc_beta_within_site <- function(df, effort, seed) {
  if (!is.null(seed)) set.seed(seed)
  
  df %>%
    group_by(YEAR, LOCALITY_ID) %>%         # iterate over each site × year
    group_split() %>%
    map_dfr(function(site_dat) {
      year_val <- site_dat$YEAR[1]
      site_val <- site_dat$LOCALITY_ID[1]
      
      # Skip sites with too few checklists
      n_checklists <- n_distinct(site_dat$SAMPLING_EVENT_IDENTIFIER)
      if (n_checklists < effort) return(NULL)
      
      # Randomly select 'effort' checklists
      chosen_chk <- site_dat %>%
        distinct(SAMPLING_EVENT_IDENTIFIER) %>%
        slice_sample(n = effort)
      
      sampled <- site_dat %>%
        semi_join(chosen_chk, by = "SAMPLING_EVENT_IDENTIFIER")
      
      # Build checklist-level community matrix
      comm <- sampled %>%
        group_by(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME) %>%
        summarise(abundance = sum(OBSERVATION_COUNT, na.rm = TRUE), .groups = "drop") %>%
        pivot_wider(names_from = COMMON_NAME, values_from = abundance, values_fill = 0) %>%
        as.data.frame()
      
      rownames(comm) <- comm$SAMPLING_EVENT_IDENTIFIER
      comm$SAMPLING_EVENT_IDENTIFIER <- NULL  # rows = checklists, cols = species
      
      # Calculate beta within this site
      beta_out <- calc_comm_div(
        comm,
        index = c("S", "S_PIE", "S_C"),
        extrapolate = TRUE,
        scales = "beta",
        C_target_gamma = 0.75
      )
      
      tibble(
        YEAR = year_val,
        LOCALITY_ID = site_val,
        scale = beta_out$scale,
        index = beta_out$index,
        sample_size = beta_out$sample_size,
        effort = beta_out$effort,
        gamma_coverage = beta_out$gamma_coverage,
        value = beta_out$value
      )
    })
}

# SUMMER
betas_within_site_summer <- calc_beta_within_site(wet_all_summer, effort = 5, seed = 5)

betas_within_site_summer <- betas_within_site_summer %>%
  mutate(SEASON = "Summer")

betas_within_site_summer <- betas_within_site_summer %>%
  filter(!is.infinite(value))

# FALL
betas_within_site_fall <- calc_beta_within_site(wet_all_fall, effort = 5, seed = 5)

betas_within_site_fall <- betas_within_site_fall %>%
  mutate(SEASON = "Fall")

betas_within_site_fall <- betas_within_site_fall %>%
  filter(!is.infinite(value))

# SPRING
betas_within_site_spring <- calc_beta_within_site(wet_all_spring, effort = 5, seed = 5)

betas_within_site_spring <- betas_within_site_spring %>%
  mutate(SEASON = "Spring")

betas_within_site_spring <- betas_within_site_spring %>%
  filter(!is.infinite(value))

# WINTER
betas_within_site_winter <- calc_beta_within_site(wet_all_winter, effort = 5, seed = 5)

betas_within_site_winter <- betas_within_site_winter %>%
  mutate(SEASON = "Winter")

betas_within_site_winter <- betas_within_site_winter %>%
  filter(!is.infinite(value))

betas_within_site <- rbind(betas_within_site_winter, betas_within_site_spring, betas_within_site_summer, betas_within_site_fall)

# saveRDS(betas_within_site, "Intermediate_data/beta_results_ecoregion_dom_removed.RDS")

# calculating the means and error bars for plotting
wet_div_error <- betas_within_site %>%
  group_by(index, SEASON) %>%
  summarize(
    mean = mean(value, na.rm = TRUE),
    lower = quantile(value, 0.025, na.rm = TRUE),
    upper = quantile(value, 0.975, na.rm = TRUE)
  )


# Plotting the results
ggplot() +
  geom_jitter(
    data = betas_within_site,
    aes(
      x = factor(index, levels = c("beta_S", "beta_S_PIE", "beta_S_C")),
      y = value,
      color = SEASON
    ),
    alpha = 0.05,
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.5)
  ) +
  geom_point(
    data = wet_div_error,
    aes(
      x = factor(index, levels = c("beta_S", "beta_S_PIE", "beta_S_C")),
      y = mean,
      color = SEASON
    ),
    size = 3,
    shape = 21,
    fill = "white",
    stroke = 1.2,
    position = position_dodge(width = 0.5)
  ) +
  geom_errorbar(
    data = wet_div_error,
    aes(
      x = factor(index, levels = c("beta_S", "beta_S_PIE", "beta_S_C")),
      ymin = lower,
      ymax = upper,
      color = SEASON
    ),
    linewidth = 1.2,
    width = 0.2,
    position = position_dodge(width = 0.5)
  ) +
  geom_hline(yintercept = 1, color = "dodgerblue", linetype = "dashed", linewidth = 1) +
  theme_bw() +
  ylab("Value") +
  xlab("Diversity Index") +
  scale_x_discrete(labels = c(
    "beta_S" = "βS",
    "beta_S_PIE" = "βSPIE",
    "beta_S_C" = "βC"
  )) +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "right")

# need to add the GIW area per state as an explanatory variable

###### building a model to explain variation in aggregation ############

# first I need to append the beta stats into a table with the env data
data_mod <- env %>%
  distinct(LOCALITY_ID, .keep_all = TRUE) %>%
  full_join(betas_within_site, by = "LOCALITY_ID")

data_mod_clean <- data_mod %>%
  filter(!is.na(value))

data_mod_clean$value <-  as.numeric(data_mod_clean$value)
data_mod_clean$YEAR <- as.factor(data_mod_clean$YEAR)


for (idx in unique(data_mod$index)) {
  cat("\n===== Index:", idx, "=====\n")
  
  df <- data_mod_clean %>% filter(index == idx)
  
  mod <- glmmTMB(
    value ~ SEASON +
      NA_L1NAME +
      (1 | YEAR) + # random intercept for YEAR
      (1 | LOCALITY_ID), # random intercept for GIW site ID, since we could have the same site on different years     
    family = Gamma(link = "log"),
    data = df
  )
  print(summary(mod))
  print(check_model(mod))
  
}

