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

# DATA
wet_all <- readRDS("Intermediate_data/all_wet_check_dominants_removed.RDS")

env <- readRDS("Data/earth_engine_env_data/env_matrix.RDS")

## Creating a dataframe with just the LOCALITY_ID and ECOREGION
## to append to the wet_all df

ecoregion_loc <- env %>%
  select(LOCALITY_ID, NA_L1NAME)

wet_all <- wet_all %>%
  inner_join(ecoregion_loc, by = "LOCALITY_ID")

######### OBJECTIVE 3 - GIW scale ############
######### Within GIW beta diversity analysis ##########

# We just need to slightly alter the function, to calculate betas WITHIN a single wetland
# to see if we have high tunrover somply between checklists at the same wetland


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
        index = c("S", "S_n", "S_PIE", "S_C"),
        extrapolate = TRUE,
        effort = 25,
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


betas_within_site <- calc_beta_within_site(wet_all, effort = 5, seed = 5)


# calculating the means and error bars for plotting
wet_div_error <- betas_within_site %>%
  group_by(index) %>%
  summarize(
    mean = mean(value, na.rm = TRUE),
    lower = quantile(value, 0.025, na.rm = TRUE),
    upper = quantile(value, 0.975, na.rm = TRUE)
  )


# Plotting the results
ggplot() +
  geom_jitter(data = betas_within_site, aes(x = index, y = value), alpha = 0.25, width = 0.2) +
  geom_point(data = wet_div_error, aes(x = index, y = mean), color = "darkred", size = 3, alpha = 1) +
  geom_errorbar(data = wet_div_error, aes(x = index, ymin = lower, ymax = upper, width = 0.2), color = "darkred", linewidth = 1) +
  geom_hline(yintercept = 1, color = "dodgerblue", linetype = "dashed", linewidth = 1) +
  theme_bw() +
  theme(legend.position = "none") +
  ylab("Value") +
  xlab("Diversity Index") +
  scale_x_discrete(labels = c(
    "beta_S" = "βS",
    "beta_S_n" = "βSn",
    "beta_S_PIE" = "βSPIE", 
    "beta_S_C" = "βC"))


# Plotting as geom_col
ggplot(wet_div_error, aes(x = index, y = mean)) +
  geom_col(fill = "skyblue") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3, color = "red") +
  geom_hline(yintercept = 1, color = "darkred", linetype = "dashed", linewidth = 1) +
  theme_bw() +
  ylab("Value") +
  xlab("Beta Diversity Index") +
  scale_x_discrete(labels = c(
    "beta_S" = "βS",
    "beta_S_n" = "βSn",
    "beta_S_PIE" = "βSPIE", 
    "beta_S_C" = "βC")) +
  theme(legend.position = "none",
        axis.title.x = element_text(
          margin = margin(t = 15)),
        axis.title.y = element_text(
          margin = margin(r = 15)))

ggsave("FigX_betas_WETLAND_scale.png", width = 6, height = 4,
       bg = "transparent")

# need to add the GIW area per state as an explanatory variable

###### building a model to explain variation in aggregation ############

# first I need to append the beta stats into a table with the env data
data_mod <- env %>%
  distinct(LOCALITY_ID, .keep_all = TRUE) %>%
  full_join(betas_within_site, by = "LOCALITY_ID")


for (idx in unique(data_mod$index)) {
  cat("\n===== Index:", idx, "=====\n")
  
  df <- data_mod %>% filter(index == idx)
  
  mod <- glm(value ~ rescale(built_wet) + rescale(water_wet) + rescale(log10(area_sqkm)) + 
               rescale(shan_wet) + rescale(evi_mean) + rescale(water_25km) + rescale(built_25km) +
               rescale(shan_gamma_25), family = poisson, data = df)
  print(summary(mod))
  print(check_model(mod))
  
}

hist(betas_within_site$value)
# wetland area
data_mod %>%
  filter(index == "beta_S") %>%
  ggplot() +
  geom_point(aes(x = log(area_sqkm*1000), y = value)) +
  geom_smooth(aes(x = log(area_sqkm*1000), y = value), method = "lm") +
  theme_bw()

# EVI
data_mod %>%
  filter(index == "beta_S") %>%
  ggplot() +
  geom_point(aes(x = evi_mean, y = value)) +
  geom_smooth(aes(x = evi_mean, y = value), method = "lm") +
  theme_bw()

# habitat heterogeneity
ggplot(data = df) +
  geom_point(aes(x = shan_wet, y = value)) +
  geom_smooth(aes(x = shan_wet, y = value), method = "lm") +
  theme_bw()

# Water 25km
data_mod %>%
  filter(index == "beta_S") %>%
  ggplot() +
  geom_point(aes(x = water_25km, y = value)) +
  geom_smooth(aes(x = water_25km, y = value), method = "lm") +
  theme_bw()
