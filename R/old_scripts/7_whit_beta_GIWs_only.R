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

# let's subset to make this managebale

env_filt <- env %>%
  select(LOCALITY_ID, NA_L1NAME)

# Combine wet_coords and env_filt
env_filt <- env_filt %>%
  inner_join(wet_coords, by = "LOCALITY_ID")

check_samples <- wet_dat %>%
  select(LOCALITY_ID, SAMPLING_EVENT_IDENTIFIER) %>%
  distinct(LOCALITY_ID, SAMPLING_EVENT_IDENTIFIER) %>%
  group_by(LOCALITY_ID) %>%
  slice_sample(n = 73)

# Now we can filter our wet_dat to only include SAMPLING_EVENT_IDENTIFIERS in
# check_smaples

###### This is the 73 checklist minimum dataset #######
wet_dat_73 <- wet_dat %>%
  filter(SAMPLING_EVENT_IDENTIFIER %in% check_samples$SAMPLING_EVENT_IDENTIFIER) %>%
  group_by(LOCALITY_ID, COMMON_NAME) %>%
  summarize(
    count = sum(OBSERVATION_COUNT)
  ) %>%
  pivot_wider(names_from = COMMON_NAME, values_from = count, values_fill = 0)

#### Now we can append the Ecoregion to this data frame to calculate gamma
#### as species richness within the ecoregion

wet_eco <- wet_dat_73 %>%
  cbind(env_filt$NA_L1NAME) %>%
  rename("ecoregion" = "...526")

wet_eco_filt <- wet_eco[ ,2:526]

# first select all the columns that are species names
species_cols <- setdiff(names(wet_eco_filt), "ecoregion")  # all species columns

# Create a species richness column, then summarize by ecoregion
eco_gamma <- wet_eco_filt %>%
  rowwise() %>%
  mutate(species_richness = sum(c_across(all_of(species_cols)) > 0)) %>%
  ungroup() %>%
  group_by(ecoregion) %>%
  summarise(mean_site_sr = mean(species_richness),
            min_sr = min(species_richness),
            gamma_sr = max(species_richness),
            n_sites = n())


##### Create the same type of dataframe by LOCALITY_ID
loc_alpha <- wet_dat_73 %>%
  rowwise() %>%
  mutate(species_richness = sum(c_across(all_of(species_cols)) > 0)) %>%
  ungroup() %>%
  group_by(LOCALITY_ID) %>%
  summarise(alpha_sr = max(species_richness))

loc_alpha <- loc_alpha %>%
  cbind(env_filt$NA_L1NAME) %>%
  rename("ecoregion" = "env_filt$NA_L1NAME")

## Now we can join the two and calculate Whitaker's Beta as gamma / alpha

beta_dat <- full_join(loc_alpha, eco_gamma, by = "ecoregion") %>%
  mutate(whit_beta_eco = gamma_sr / alpha_sr,
         tot_gamma = 524,
         whit_beta_tot = 524 / alpha_sr)


## Now we can make a plot of the whitaker's beta values

ggplot() +
  geom_jitter(data = beta_dat, aes(x = ecoregion, y = whit_beta_eco, color = ecoregion)) +
  theme_bw() +
  xlab("Ecoregion") +
  ylab("Whitaker's Beta") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none")
