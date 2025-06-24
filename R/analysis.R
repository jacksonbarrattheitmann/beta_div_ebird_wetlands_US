# analysis script

library(dplyr)
library(tidyverse)
library(vegan)
library(mobr)
library(ggplot2)
library(betapart)
library(lme4)
library(performance)

######### Within GIW beta diversity analysis ##########

#read in env and wet_comm data
env <- readRDS("Intermediate_data/env_all_wet_check.RDS")

wet_dat <- readRDS("Intermediate_data/all_wet_check_summarized.RDS")

## make SAMPLING_EVENT_IDENTIFIER into the row names

env <- column_to_rownames(env, var = "SAMPLING_EVENT_IDENTIFIER")

# do the same for wet_dat
wet_dat <- column_to_rownames(wet_dat, var = "SAMPLING_EVENT_IDENTIFIER")

# need to remove LOCALITY_ID column
wet_dat <- wet_dat[, 2:529]

######### now time to run mobr ##################
calc_C_target(wet_dat)

wet_div <- tibble(wet_dat) %>% 
  group_by(group = env$LOCALITY_ID) %>% 
  group_modify(~ calc_comm_div(.x, index = c('S','S_C'),
                               extrapolate = TRUE), scales = c("beta"), C_target_gamma = 0.75)

# the defualt mobr plot
plot_comm_div(wet_div)

# subsetting just the beta_S_C values for plotting in ggplot
wet_div_betas <- wet_div %>%
  filter(index == "beta_S_C")
# creating a table for the errorbar values
wet_div_error <- wet_div %>%
  filter(index == "beta_S_C")  %>%
  select(value, index) %>%
  group_by(index) %>%
  summarize(median = median(value),
            lower = quantile(value, probs = 0.025),
            upper = quantile(value, probs = 0.975))

# Plotting the results
ggplot() +
  geom_jitter(data = wet_div_betas, aes(x = index, y = value), alpha = 0.5, width = 0.2) +
  geom_point(data = wet_div_betas, aes(x = index, y = median(value), color = "red")) +
  geom_errorbar(data = wet_div_error, aes(x = index, ymin = lower, ymax = upper, width = 0.1), color = "red") +
  theme_bw() +
  theme(legend.position = "none") +
  ylab("Value")


###### building a model to explain variation in aggregation ############
wet_div_betas_filt <- wet_div_betas %>%
  rename(LOCALITY_ID = group)

# first I need to append the beta stats into a table with the env data
data_mod <- env %>%
  distinct(LOCALITY_ID, .keep_all = TRUE) %>%
  full_join(wet_div_betas_filt, by = "LOCALITY_ID")


mod <- glm(value ~ built_wet + flooded_vegetation_wet + area_sqkm + shan_wet + trees_wet, 
           family = quasipoisson(), data = data_mod)
summary(mod)
check_model(mod)

plot(log(data_mod$area_sqkm), data_mod$value)
