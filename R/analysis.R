# analysis script

library(dplyr)
library(tidyverse)
library(vegan)
library(mobr)
library(ggplot2)
library(betapart)

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
calc_C

wet_div <- tibble(wet_dat) %>% 
  group_by(group = env_filt$LOCALITY_ID) %>% 
  group_modify(~ calc_comm_div(.x, index = c('N','S','S_n', 'S_PIE', 'S_C'), effort = 25,
                               extrapolate = TRUE), scales = c("alpha", "gamma", "beta"))
