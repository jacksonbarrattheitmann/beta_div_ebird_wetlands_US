# Creating the species-site matrix 
# this script is the script to calculate wetland specific measures of diversity
# both richness and abundance
# in 'total', defined as: XXXXXXXXXXX
# in mean across checklists, defined as: XXXXXXXXXXx
# and mean/median resampled with 50 checklists at each wetland, defined as: XXXXXXX

library(dplyr)
library(tidyverse)


file_list <- list.files(path = "Data/eBird_local_alpha_level", pattern = "\\.RDS", full.names = TRUE)
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

