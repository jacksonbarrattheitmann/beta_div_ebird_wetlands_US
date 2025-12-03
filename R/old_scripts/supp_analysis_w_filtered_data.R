### Analysis script with just 
### checklists from the summer of 2021 that were stationary, and less
### than 60 minutes

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
library(geosphere)
library(purrr)
library(sf)
library(units)
library(multcompView)
library(broom)
########## WETLAND scale ##########

## first let's read in the checklist level eBird data that's filtered
wet_dat <- readRDS("Intermediate_data/all_wet_check_filt_2021_summer.RDS")

#for simplicity let's just filter down to only 5 checklists for every site
check_samples <- wet_dat %>%
  select(LOCALITY_ID, SAMPLING_EVENT_IDENTIFIER) %>%
  distinct(LOCALITY_ID, SAMPLING_EVENT_IDENTIFIER) %>%
  group_by(LOCALITY_ID) %>%
  slice_sample(n = 5)

# Need an intermediate df with just the sampling event identifiers and LOCALITY_ID
wet_check <- check_samples %>%
  select(LOCALITY_ID, SAMPLING_EVENT_IDENTIFIER)

# now read in the env dataframe
env <- readRDS("Data/earth_engine_env_data/env_matrix.RDS")

# now append
wet_check <- wet_check %>%
  left_join(env, by = "LOCALITY_ID") %>%
  distinct(SAMPLING_EVENT_IDENTIFIER, .keep_all = TRUE)

env_mat <- column_to_rownames(wet_check, var = "SAMPLING_EVENT_IDENTIFIER")

# Make this into a site-species matrix, except now rows
# are checklists instead of sites

wet_dat_mat <- wet_dat %>%
  filter(SAMPLING_EVENT_IDENTIFIER %in% check_samples$SAMPLING_EVENT_IDENTIFIER) %>%
  group_by(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME) %>%
  summarize(
    count = sum(OBSERVATION_COUNT)
  ) %>%
  pivot_wider(names_from = COMMON_NAME, values_from = count, values_fill = 0)


# do the same for wet_dat
wet_dat_mat <- column_to_rownames(wet_dat_mat, var = "SAMPLING_EVENT_IDENTIFIER")



######### now time to run mobr ##################

wet_div <- tibble(wet_dat_mat) %>%
  group_by(group = env_mat$LOCALITY_ID) %>%
  group_modify(~ calc_comm_div(.x,
                               index = c('S', 'S_n','S_PIE', 'S_C'),
                               extrapolate = TRUE,
                               effort = 25,
                               scales = "beta",
                               C_target_gamma = 0.75))


wet_div_error <- wet_div %>%
  group_by(index) %>%
  summarize(
    mean = mean(value, na.rm = TRUE),
    lower = quantile(value, 0.025, na.rm = TRUE),
    upper = quantile(value, 0.975, na.rm = TRUE)
  )


# Plotting the results
ggplot() +
  geom_jitter(data = wet_div, aes(x = index, y = value), alpha = 0.25, width = 0.2) +
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


####### CONTINETNAL ########

## Now we can also do this at the CONTINENTAL scale

# So let's collapse and aggregate the data by LOCALITY_ID 
# instead of by SAMPLING_EVENT_IDENTIFIER

wet_dat_cont <- wet_dat %>%
  filter(SAMPLING_EVENT_IDENTIFIER %in% check_samples$SAMPLING_EVENT_IDENTIFIER) %>%
  group_by(LOCALITY_ID, COMMON_NAME) %>%
  summarize(
    count = sum(OBSERVATION_COUNT)
  ) %>%
  pivot_wider(names_from = COMMON_NAME, values_from = count, values_fill = 0)

# Need an intermediate df with just the LOCALITY_ID
wet_check_cont <- wet_dat_cont %>%
  select(LOCALITY_ID)


# now append
wet_check_cont <- wet_check_cont %>%
  left_join(env, by = "LOCALITY_ID")

env_mat_cont <- column_to_rownames(wet_check_cont, var = "LOCALITY_ID")

# do the same for wet_dat
wet_dat_cont <- column_to_rownames(wet_dat_cont, var = "LOCALITY_ID")


## The mobR approach
betas_ALL <- calc_comm_div(wet_dat_cont, index = c('S','S_n','S_PIE', 'S_C'),
                           extrapolate = TRUE, effort = 25, scales = c("beta"), C_target_gamma = 0.75)

########## PLOTS for OBJECTIVE 1 ##################

ggplot(data = betas_ALL) +
  geom_col(aes(x = index, y = value), fill = "skyblue") + theme_bw() +
  geom_hline(yintercept = 1, color = "darkred", linetype = "dashed", linewidth = 1) +
  xlab("Beta Diversity Index") +
  ylab("Value") +
  scale_x_discrete(  labels = c(
    "beta_S" = "βS",
    "beta_S_n" = "βSn",
    "beta_S_PIE" = "βSPIE", 
    "beta_S_C" = "βC"
  )) +
  theme(legend.position = "none",
        axis.title.x = element_text(
          margin = margin(t = 15)),
        axis.title.y = element_text(
          margin = margin(r = 15)))

######### ECOREGION ########

env_mat_cont <- rownames_to_column(env_mat_cont, var = "LOCALITY_ID")

# need spatial coords for the LOCALITY_IDs
wet_coords <- readRDS("Intermediate_data/locality_ids_long_lat.RDS")

# should join the coords with the env dataframe
env_mat_cont <- env_mat_cont %>%
  inner_join(wet_coords, by = "LOCALITY_ID")

beta_results <- calculate_beta_ecoregion_sf(wet_dat_cont, env_mat_cont, n_reps = 50)

## only have EASTERN TEMPERATE FORESTS due to sampling issues

beta_div_error <- beta_results %>%
  group_by(index) %>%
  summarize(
    mean = mean(value, na.rm = TRUE),
    lower = quantile(value, 0.025, na.rm = TRUE),
    upper = quantile(value, 0.975, na.rm = TRUE)
  )


# Plotting the results
ggplot() +
  geom_jitter(data = beta_results, aes(x = index, y = value), alpha = 0.25, width = 0.2) +
  geom_point(data = beta_div_error, aes(x = index, y = mean), color = "darkred", size = 3, alpha = 1) +
  geom_errorbar(data = beta_div_error, aes(x = index, ymin = lower, ymax = upper, width = 0.2), color = "darkred", linewidth = 1) +
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
