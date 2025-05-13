# Creating the species-site matrix 
# this script is the script to calculate wetland specific measures of diversity
# both richness and abundance
# in 'total', defined as: XXXXXXXXXXX
# in mean across checklists, defined as: XXXXXXXXXXx
# and mean/median resampled with 50 checklists at each wetland, defined as: XXXXXXX

library(dplyr)
library(tidyverse)

file_list <- list.files(path = "Data/eBird_local_alpha_level/", pattern = "\\.RDS", full.names = TRUE)
r <- lapply(file_list, readRDS)
s <- bind_rows(r, .id = "column_label")
s <- as.data.frame(s)

s$OBSERVATION_COUNT <- as.numeric(s$OBSERVATION_COUNT)
s$COMMON_NAME <- as.factor(s$COMMON_NAME)
s$SAMPLING_EVENT_IDENTIFIER <- as.factor(s$SAMPLING_EVENT_IDENTIFIER)

#test <- s %>% 
#  filter(SAMPLING_EVENT_IDENTIFIER == "S13774755")

#grabbing data for Muscovy Duck and Rock Piegon
# because of weird issue Jackson add details
temp <- s %>% dplyr::filter(CATEGORY == "domestic") %>%
  filter(COMMON_NAME == "Rock Pigeon" | COMMON_NAME == "Muscovy Duck", OBSERVATION_COUNT >= 1)

# Combining the rock pigeon and muscovy duck df with the larger s df
s_joined <- s %>%
  dplyr::filter(CATEGORY %in% c("species", "issf"), OBSERVATION_COUNT >= 1) %>%
  bind_rows(temp)

get_group_checklists <- s_joined %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, GROUP_IDENTIFIER) %>%
  distinct() %>%
  replace_na(list(GROUP_IDENTIFIER="not shared")) %>%
  dplyr::filter(GROUP_IDENTIFIER != "not shared") %>%
  group_by(GROUP_IDENTIFIER) %>%
  sample_n(1)

s_not_shared_checklists <- s_joined %>%
  replace_na(list(GROUP_IDENTIFIER="not shared")) %>%
  dplyr::filter(GROUP_IDENTIFIER == "not shared")

s_trimmed_group_shared_checklists <- s_joined %>%
  dplyr::filter(SAMPLING_EVENT_IDENTIFIER %in% local(get_group_checklists$SAMPLING_EVENT_IDENTIFIER))

s_all_data_cleaned <- s_not_shared_checklists %>%
  bind_rows(s_trimmed_group_shared_checklists)

# filtering by species and issf, and Observation greater than or equal to 1 to get rid of Xs
species_counts_fixed <- s_all_data_cleaned %>%
  group_by(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME) %>%
  summarize(OBSERVATION_COUNT=sum(OBSERVATION_COUNT))

lists_meta <- s_all_data_cleaned %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, 6:22) %>%
  dplyr::select(-OBSERVER_ID, -GROUP_IDENTIFIER, -CATEGORY) %>%
  distinct()

s_final <- lists_meta %>%
  left_join(species_counts_fixed) 

# WRITING RDS file with all the wetland data raw, before summarizing

saveRDS(s_final, "Data/eBird_local_alpha_level/ebird_alpha_wetlands_raw.RDS")


# Calculating total richness, total checklists, total observers, total abundance
s_matrix <- s_final %>%
  group_by(LOCALITY_ID) %>%
  summarise(total_checklists = length(unique(SAMPLING_EVENT_IDENTIFIER)), 
            total_richness = length(unique(COMMON_NAME)), 
            total_observers = sum(NUMBER_OBSERVERS),
            total_N = sum(OBSERVATION_COUNT),
            group = "All Species")

# calculating mean richness and abundance
# across all checklists submitted at a locality_id
mean_diversity_across_checklists <- s %>%
  dplyr::filter(CATEGORY == "species") %>%
  group_by(LOCALITY_ID, SAMPLING_EVENT_IDENTIFIER) %>%
  summarize(checklist_rich=length(unique(COMMON_NAME)),
            checklist_abund=sum(OBSERVATION_COUNT)) %>%
  ungroup() %>%
  group_by(LOCALITY_ID) %>%
  summarize(mean_richness_across_checklists=mean(checklist_rich),
            median_richness_across_checklists=median(checklist_rich),
            mean_abund_across_checklists=mean(checklist_abund),
            median_abund_across_checklists=median(checklist_abund))


# calculate a resampled species richness
# downsampling to 50 checklists at each wetland
downsample_function <- function(wetland_id, number_of_checklists){
  
  message(paste0("calculating richness for ", wetland_id))
  
  resample_function <- function(draw_number){
    
    wetland_dat <- s %>%
      dplyr::filter(LOCALITY_ID == wetland_id)
    
    checklist_sample <- wetland_dat %>%
      dplyr::select(SAMPLING_EVENT_IDENTIFIER) %>%
      distinct() %>%
      sample_n(number_of_checklists)
    
    sample_richness <- wetland_dat %>%
      dplyr::filter(SAMPLING_EVENT_IDENTIFIER %in% checklist_sample$SAMPLING_EVENT_IDENTIFIER) %>%
      group_by(LOCALITY_ID) %>%
      summarize(sample_richness=length(unique(COMMON_NAME)),
                sample_abund=sum(OBSERVATION_COUNT)) %>%
      mutate(draw=draw_number)
    
    return(sample_richness)
  }
  
  resampled_richness <- bind_rows(lapply(1:100, function(x) resample_function(x)))
                                  
  final_locality_resampled_rich <- resampled_richness %>%
    group_by(LOCALITY_ID) %>%
    summarize(resampled_median_richness=median(sample_richness),
              resampled_mean_richness=mean(sample_richness),
              resampled_median_abund=median(sample_abund),
              resampled_mean_abund=mean(sample_abund))
  
  return(final_locality_resampled_rich)
}

resampled_richness_per_locality <- bind_rows(lapply(unique(s$LOCALITY_ID), function(x) downsample_function(x, 50)))

# appending this output to the s_matrix data frame
s_final <- s_matrix %>%
  inner_join(mean_diversity_across_checklists, by="LOCALITY_ID") %>%
  inner_join(resampled_richness_per_locality, by = "LOCALITY_ID")


# Still need to work on the resample for the median richness, from 100 checklists
saveRDS(s_final, file = "Data/species_matrix/s_matrix.RDS")


############# Some summary statistics


# total number of observations
nrow(s_final)

# total number of birds observed
sum(s_final$OBSERVATION_COUNT)

# total number of unqiue species
length(unique(s_final$COMMON_NAME))

# total number of observers
length(unique(s_final$SAMPLING_EVENT_IDENTIFIER))

# average number of observers per wetland
length(unique(s_final$SAMPLING_EVENT_IDENTIFIER))/207

#number of checklists at fewest and most vistied wetlands
w_most <- s_final %>%
  filter(LOCALITY_ID == "L208918") %>%
  count(SAMPLING_EVENT_IDENTIFIER)

w_least <- s_final %>%
  filter(LOCALITY_ID == "L11338314") %>%
  count(SAMPLING_EVENT_IDENTIFIER)

