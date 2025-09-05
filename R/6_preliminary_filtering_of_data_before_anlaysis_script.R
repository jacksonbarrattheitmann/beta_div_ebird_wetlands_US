# Generating single site-species matrix for all 207 wetlands

library(dplyr)
library(tidyverse)
library(vegan)
library(mobr)
library(ggplot2)
library(betapart)
library(lubridate)
library(ggrepel)

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

# Max = 22650
# Min = 92

###### When calculating beta diversity we're going to have to deal with the issue
###### of temporal beta when calculating spatial beta
###### so let's just summarize the data into a checklist x species matrix
###### so that we can sub-sample/bootstrap only within years and seasons

### So let's take wet_dat
### and save it first just by season, across all years, that way we can
### write a function to deal with year in the analysis script

# first let's split up OBSERVATION_DATE into 3 separate columns for easy filtering
wet_dat <- wet_dat %>%
  mutate(
    OBSERVATION_DATE = ymd(OBSERVATION_DATE),
    YEAR  = year(OBSERVATION_DATE),
    MONTH = month(OBSERVATION_DATE),
    DAY   = day(OBSERVATION_DATE)
  )

wet_dat <- wet_dat %>%
  mutate(
    SEASON = case_when(
      MONTH %in% c(12, 1, 2) ~ "WINTER",
      MONTH %in% c(3, 4, 5)  ~ "SPRING",
      MONTH %in% c(6, 7, 8)  ~ "SUMMER",
      MONTH %in% c(9, 10, 11) ~ "FALL"
    )
  )

## Now let's see if we can't create a filter to have at least 5 checklists
## in only the summer months, in any given year

wet_dat_analysis <- wet_dat %>%
  filter(SEASON == "SUMMER") %>%
  group_by(LOCALITY_ID, YEAR) %>%
  mutate(num_check = length(unique(SAMPLING_EVENT_IDENTIFIER)))

## here's the 5 checklist filter
wet_dat_analysis <- wet_dat_analysis %>%
  filter(num_check >= 5)


## Now let's look at how many LOCALITY_IDs we'll have
length(unique(wet_dat_analysis$LOCALITY_ID))

# 195 with at least 5 checklists during one summer season in at least 1 year
# let's save this as an RDS for our analysis script

## Let's make sure we have independence in the sampling units
# by grabbing only a single checklist from a single LOCALITY_ID on a single DAY
wet_dat_poss <- wet_dat_analysis %>%
  distinct(LOCALITY_ID, SAMPLING_EVENT_IDENTIFIER, OBSERVATION_DATE)

# looks like this is not a problem so we can go ahead with our current 
# dataset in wet_dat_anlaysis
wet_dat_analysis_2 <- wet_dat_analysis %>%
  filter(SAMPLING_EVENT_IDENTIFIER %in% wet_dat_poss$SAMPLING_EVENT_IDENTIFIER)

## SAVING THE FILE
#saveRDS(wet_dat_analysis, "Intermediate_data/all_wet_check_for_analysis.RDS")


#### Summary statistics
spp_dat <- wet_dat_analysis %>%
  group_by(COMMON_NAME) %>%
  summarise(abundance = sum(OBSERVATION_COUNT, na.rm = TRUE))

# plot for supplement
ggplot(data = spp_dat) +
  geom_histogram(aes(x = abundance), bins = 30, fill = "lightblue") +
  xlab("Total Abundance") +
  ylab("Count of species") +
  scale_x_continuous(labels = label_comma()) +
  theme_bw()

## total species

length(unique(spp_dat$COMMON_NAME))
#  497


## total abundance

sum(spp_dat$abundance)
# 4,053,763

####
# 162592

# Percentage of sites with the species present
spp_dat_2 <- wet_dat_analysis %>%
  distinct(LOCALITY_ID, COMMON_NAME) %>%
  count(COMMON_NAME, name = "n_sites") %>%
  mutate(total_sites = n_distinct(wet_dat_analysis$LOCALITY_ID),
         prop_sites = n_sites / total_sites)

## Plot for supplement

spp_dat <- spp_dat %>%
  inner_join(spp_dat_2, by = "COMMON_NAME")

ggplot(data = spp_dat, aes(x = prop_sites, y = abundance)) +
  geom_point() +
  geom_text_repel(aes(label = COMMON_NAME)) +
  scale_y_continuous(labels = label_comma()) +
  xlab("Proportion of sites with species present") +
  ylab("Total Abundance of species") +
  theme_bw()



######### SUPPLEMENTRAY DATASET ############
##### Let's also create a secondary dataset for analysis that excludes
##### our 5 most abundant species, to see if rarity is better explained by
##### env variation

spp_dat_filt <- spp_dat %>%
  filter(abundance < 162592)

### Now filter the wet_dat_analysis dataframe for only these species

wet_dat_no_dom <- wet_dat_analysis %>%
  filter(COMMON_NAME %in% spp_dat_filt$COMMON_NAME)

#saveRDS(wet_dat_no_dom, "Intermediate_data/all_wet_check_dominants_removed.RDS")

#### Supplementary Analysis with just one year, and only stationary checklists
# let's deal with temporal bias + distance traveled
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

#saveRDS(sum_check_stat_filt, "Intermediate_data/all_wet_check_filt_2021_summer.RDS")



