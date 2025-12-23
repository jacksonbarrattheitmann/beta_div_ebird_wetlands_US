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

######## SUMMER ###########
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
saveRDS(wet_dat_analysis, "Intermediate_data/all_wet_check_for_analysis_SUMMER.RDS")



###### WINTER #######

wet_dat_analysis <- wet_dat %>%
  filter(SEASON == "WINTER") %>%
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
#saveRDS(wet_dat_analysis, "Intermediate_data/all_wet_check_for_analysis_WINTER.RDS")


###### SPRING #########
wet_dat_analysis <- wet_dat %>%
  filter(SEASON == "SPRING") %>%
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
#saveRDS(wet_dat_analysis, "Intermediate_data/all_wet_check_for_analysis_SPRING.RDS")


###### FALL ######

wet_dat_analysis <- wet_dat %>%
  filter(SEASON == "FALL") %>%
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
#saveRDS(wet_dat_analysis, "Intermediate_data/all_wet_check_for_analysis_FALL.RDS")


## LOAD in the data files to remove hyperabundants, and make 
## species plots by season

wet_dat_sum <- readRDS("Intermediate_data/all_wet_check_for_analysis_SUMMER.RDS")
wet_dat_win <- readRDS("Intermediate_data/all_wet_check_for_analysis_WINTER.RDS")
wet_dat_fall <- readRDS("Intermediate_data/all_wet_check_for_analysis_FALL.RDS")
wet_dat_spr <- readRDS("Intermediate_data/all_wet_check_for_analysis_SPRING.RDS")

wet_dat_all <- rbind(wet_dat_spr, wet_dat_win, wet_dat_sum, wet_dat_fall)

#total checklists
length(unique(wet_dat_sum$SAMPLING_EVENT_IDENTIFIER))

# total sites
length(unique(wet_dat_sum$LOCALITY_ID))

#total species
length(unique(wet_dat_all$COMMON_NAME))

#total birds counted
sum(wet_dat_all$OBSERVATION_COUNT)

#### Summary statistics
spp_dat <- wet_dat_all %>%
  group_by(COMMON_NAME, SEASON) %>%
  summarise(abundance = sum(OBSERVATION_COUNT, na.rm = TRUE))

# total sites by season
spp_dat_season <- wet_dat_all %>%
  group_by(SEASON) %>%
  summarize(total_sites = n_distinct(LOCALITY_ID))

median(spp_dat$abundance)
min(spp_dat$abundance)
max(spp_dat$abundance)


# checklists level data


# Percentage of sites with the species present
spp_dat_2 <- wet_dat_all %>%
  group_by(SEASON) %>%
  distinct(LOCALITY_ID, COMMON_NAME) %>%
  count(COMMON_NAME, name = "n_sites") %>%
  inner_join(spp_dat_season, by = "SEASON") %>%
  mutate(prop_sites = n_sites / total_sites)

## Plot for supplement

spp_dat <- spp_dat %>%
  inner_join(spp_dat_2, by = c("SEASON", "COMMON_NAME"))

ggplot(data = spp_dat, aes(x = prop_sites, y = abundance, color = SEASON)) +
  geom_point() +
  geom_text_repel(aes(label = COMMON_NAME)) +
  scale_y_continuous(labels = label_comma()) +
  facet_wrap(~SEASON) +
  xlab("Proportion of sites with species present") +
  ylab("Total Abundance of species") +
  theme_bw() +
  scale_color_brewer(palette = "Set1")



######### SUPPLEMENTRAY DATASET ############
##### Let's also create a secondary data set for analysis that excludes
##### our 5 most abundant species, to see if rarity is better explained by
##### env variation

# Let's find the 5 most abundant species
# within each season, and remove those

most_abun <- spp_dat %>%
  group_by(SEASON) %>%
  arrange(desc(abundance)) %>%
  slice_head(n = 5)

### Now filter the wet_dat_analysis dataframe for only these species

most_abun_sum <- most_abun %>%
  filter(SEASON == "SUMMER")

most_abun_win <- most_abun %>%
  filter(SEASON == "WINTER")

most_abun_spr <- most_abun %>%
  filter(SEASON == "SPRING")

most_abun_fall <- most_abun %>%
  filter(SEASON == "FALL")

## No we can filter each of the season dfs

wet_dat_sum_dom_rem <- wet_dat_sum %>%
  filter(!COMMON_NAME %in% most_abun_sum$COMMON_NAME)

saveRDS(wet_dat_sum_dom_rem, "Intermediate_data/all_wet_check_dominants_removed_summer.RDS")

wet_dat_win_dom_rem <- wet_dat_win %>%
  filter(!COMMON_NAME %in% most_abun_win$COMMON_NAME)

saveRDS(wet_dat_win_dom_rem, "Intermediate_data/all_wet_check_dominants_removed_winter.RDS")

wet_dat_spr_dom_rem <- wet_dat_spr %>%
  filter(!COMMON_NAME %in% most_abun_spr$COMMON_NAME)

saveRDS(wet_dat_spr_dom_rem, "Intermediate_data/all_wet_check_dominants_removed_spring.RDS")

wet_dat_fall_dom_rem <- wet_dat_fall %>%
  filter(!COMMON_NAME %in% most_abun_fall$COMMON_NAME)

saveRDS(wet_dat_fall_dom_rem, "Intermediate_data/all_wet_check_dominants_removed_fall.RDS")








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



