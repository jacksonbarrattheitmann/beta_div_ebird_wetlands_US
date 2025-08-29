# Generating single site-species matrix for all 207 wetlands

library(dplyr)
library(tidyverse)
library(vegan)
library(mobr)
library(ggplot2)
library(betapart)

# First let's read in all the GIW data as a single dataframe

wet_dat <- readRDS("Data/eBird_local_alpha_level/ebird_alpha_wetlands_raw.RDS")

# The ENV dataset

env <- readRDS("Data/earth_engine_env_data/env_matrix.RDS")

wet_coords <- wet_dat %>%
  select(LOCALITY_ID, LONGITUDE, LATITUDE) %>%
  distinct()

#saveRDS(wet_coords, "Intermediate_data/locality_ids_long_lat.RDS")

# let's subset to make this managebale

env_filt <- env %>%
  select(LOCALITY_ID, NA_L1NAME)

# Combine wet_coords and env_filt
env_filt <- env_filt %>%
  inner_join(wet_coords, by = "LOCALITY_ID")


# let's create a site-species matrix for all 207 GIWs from the minimum number of 
# checklists that we have from each GIW

# First, let's check and see how many checklists we have at each GIW

sum_check <- wet_dat %>%
  group_by(LOCALITY_ID) %>%
  summarize(num_check = length(unique(SAMPLING_EVENT_IDENTIFIER)))

# Max = 3550
# Min = 73



# We have a few options
# 1) grab only 73 checklists from all the GIWs, and collapse them all, calculate beta 1 time
# 2) Bootstrapping approach, randomly sub-sample 25-50 and calculate beta each time

# Let's just try #1 for now, and use mobr to analyze

# First thing is we need to create a sample of SAMPLING_EVENT_IDENTIFIERS
# with length 73, for each LOCALITY_ID (eBird hotspot)

check_samples <- wet_dat %>%
  select(LOCALITY_ID, SAMPLING_EVENT_IDENTIFIER) %>%
  distinct(LOCALITY_ID, SAMPLING_EVENT_IDENTIFIER) %>%
  group_by(LOCALITY_ID) %>%
  slice_sample(n = 73)

# check to make sure this worked
# if this is an empty data frame we are good to go
test <- check_samples %>%
  group_by(LOCALITY_ID) %>%
  summarize(num = length(unique(SAMPLING_EVENT_IDENTIFIER))) %>%
  filter(num != 73)

# Here is the dataframe with all the checklist data
all_wet_check_73 <- wet_dat %>%
  filter(SAMPLING_EVENT_IDENTIFIER %in% check_samples$SAMPLING_EVENT_IDENTIFIER) %>%
  group_by(SAMPLING_EVENT_IDENTIFIER, LOCALITY_ID, COMMON_NAME) %>%
  summarize(
    count = sum(OBSERVATION_COUNT)
  ) %>%
  pivot_wider(names_from = COMMON_NAME, values_from = count, values_fill = 0)


saveRDS(all_wet_check_73, "Intermediate_data/all_wet_check_summarized.RDS")


# Need to create a env_file that had the same number of rows
env_all_wet_check <- all_wet_check_73 %>%
  select(LOCALITY_ID) %>%
  left_join(env, by = "LOCALITY_ID")

saveRDS(env_all_wet_check, "Intermediate_data/env_all_wet_check.RDS")

# Now we can filter our wet_dat to only include SAMPLING_EVENT_IDENTIFIERS in
# check_samples
wet_dat_73 <- wet_dat %>%
  filter(SAMPLING_EVENT_IDENTIFIER %in% check_samples$SAMPLING_EVENT_IDENTIFIER) %>%
  group_by(LOCALITY_ID, COMMON_NAME) %>%
  summarize(
    count = sum(OBSERVATION_COUNT)
  ) %>%
 pivot_wider(names_from = COMMON_NAME, values_from = count, values_fill = 0)

# Now we need to filter out the ecoregions with less than 10 GIWs

env_filt_10_GIWs <- env %>%
  group_by(NA_L1NAME) %>%
  filter(NA_L1NAME == "EASTERN TEMPERATE FORESTS" | NA_L1NAME == "GREAT PLAINS" |
           NA_L1NAME == "MARINE WEST COAST FOREST" | NA_L1NAME == "MEDITERRANEAN CALIFORNIA" |
           NA_L1NAME == "NORTH AMERICAN DESERTS")

wet_comm_10_GIWs <- wet_dat_73 %>%
  filter(LOCALITY_ID %in% env_filt_10_GIWs$LOCALITY_ID)

saveRDS(env_filt_10_GIWs, "Intermediate_data/env_ecoregions_10_GIWs.RDS")
saveRDS(wet_dat_73, "Intermediate_data/wet_comm_all_summarized.RDS")
saveRDS(wet_comm_10_GIWs, "Intermediate_data/wet_comm_ecoregion_10_summarized.RDS")


#### Let's subset these data few different ways for our 
#### analysis besides just using all 73 checklists aggregated

# first let's deal with temporal bias
# so let's deal with just a single summer of data 2021
## and by filtering to only stationary checklists
sum_check_stat <- wet_dat %>%
  group_by(LOCALITY_ID) %>%
  filter(PROTOCOL_TYPE == "Stationary") %>%
  filter(OBSERVATION_DATE >= as.Date("2021-05-01") & OBSERVATION_DATE <= as.Date("2021-08-30")) %>%
  filter(DURATION_MINUTES < 60) %>%
  mutate(num_check = length(unique(SAMPLING_EVENT_IDENTIFIER))) 

# Now let's remove LOCALITY_IDs that don't have at least 5 checklists
# as that will be our bootstrap size
sum_check_stat_filt <- sum_check_stat %>%
  filter(num_check >= 5)


##### Now we have a dataframe of checklists that meet our strict requirments #####
##### Now we can save it and make a function to do our analysis in mobr ######

saveRDS(sum_check_stat_filt, "Intermediate_data/all_wet_check_filt_2021_summer.RDS")



