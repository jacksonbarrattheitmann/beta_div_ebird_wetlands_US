library(dplyr)
library(tidyverse)
library(ggplot2)
library(taxize)
library(purrr)

### Checking out the Abundance data

wet_eco <- readRDS("Intermediate_data/wet_comm_ecoregion_10_summarized.RDS") %>%
  column_to_rownames(var = "LOCALITY_ID")

# the ENV data with ECOREGION = NA_L1NAME column
env_eco <- readRDS("Intermediate_data/env_ecoregions_10_GIWs.RDS") 

# AVONET table for taxonomy and functional traits
tax_func <- read.csv("Data/functional_traits_data_AVONET/AVONET1_BirdLife.csv") %>%
  rename(SCIENTIFIC_NAME = Species1)

## Just a plot of the most abundant species
species_sums <- wet_eco %>%
  summarise(across(everything(), sum, na.rm = TRUE)) %>%
  pivot_longer(cols = everything(), names_to = "species", values_to = "total") %>%
  arrange(desc(total))  # Sort from high to low

# to see the most abundant species
print(species_sums)

#
species_sums <- species_sums %>%
  mutate(pct_tot = total/sum(total))

## Let's try and get the taxonomy for all these species so we can sum them 
# different based on taxonomy

# First I need a list of the common names 
# and scientific names from the raw data

raw_dat <- readRDS("Data/eBird_local_alpha_level/ebird_alpha_wetlands_raw.RDS")
raw_dat <- raw_dat %>%
  distinct(COMMON_NAME, .keep_all = TRUE) %>%
  select(COMMON_NAME, SCIENTIFIC_NAME) %>%
  rename(species = COMMON_NAME)


## Now append this to species_sums

species_sums <- species_sums %>%
  inner_join(raw_dat, by = "species") 

species_sums <- species_sums %>%
  inner_join(tax_func, by = "SCIENTIFIC_NAME")

species_sums <- species_sums %>%
  group_by(Order1) %>%
  mutate(tot_order = sum(total),
         pct_order = sum(pct_tot))

species_sums <- species_sums %>%
  group_by(Family1) %>%
  mutate(tot_fam = sum(total),
         pct_fam = sum(pct_tot))




## Family plot of abundance
species_sums %>%
  distinct(Family1, .keep_all = TRUE) %>%
ggplot() +
  geom_col(aes(x = Family1, y = tot_fam, fill = Family1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none") +
  xlab("Family") +
  ylab("Total Abundance (N)")

## Order plot of abundance
species_sums %>%
  distinct(Order1, .keep_all = TRUE) %>%
ggplot() +
  geom_col(aes(x = Order1, y = tot_order, fill = Order1)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none") +
  xlab("Order") +
  ylab("Total Abundance (N)")





###### Doing a rank count of abundance since
###### Some a few species are super high

wet_eco <- wet_eco %>% 
  mutate(site = rownames(.))

## Taking this from wide to long format, and calculating the number
## of species occurrences
N_max <- wet_eco %>%
  pivot_longer(
    cols = -site,
    names_to = "species",
    values_to = "count"
  ) %>%
  group_by(site) %>%
  mutate(total_N = sum(count)) %>%
  slice_max(count, n = 5, with_ties = FALSE) %>%
  ungroup()

N_max <- N_max %>%
  group_by(site) %>%
  mutate(rank = rank(-count, ties.method = "first"),
         N_prop = count/total_N) %>%
  ungroup()

N_max_test <- N_max %>%
  mutate(rank_score = 6 - rank) %>% 
  # Sum scores by species
  group_by(species) %>%
  mutate(rank_sum = sum(rank_score)) %>%
  ungroup()

N_max_test <- N_max_test %>%
  rename(LOCALITY_ID = site)

## Now we can append the env df

N_max_test <- N_max_test %>%
  inner_join(env_eco, by = "LOCALITY_ID")

## Now I'm going to calculate the rank score
## again but within Ecoregion

N_max_test <- N_max_test %>%
  group_by(species, NA_L1NAME) %>%
  mutate(rank_sum_eco = sum(rank_score), 
         eco_N = sum(count)) %>%
  ungroup()

## Now keep just the distinct rows for plotting

N_filt <- N_max_test %>%
  distinct(species, NA_L1NAME, .keep_all = TRUE)


N_filt <- N_filt %>%
  mutate(species = fct_reorder(species, rank_sum_eco, .fun = max, .desc = TRUE))

### Now finally, we can make some plots

N_filt %>%
  filter(rank == 1) %>%
ggplot() +
  geom_col(aes(x = species, y = rank_sum_eco, fill = NA_L1NAME)) +
  facet_wrap(~NA_L1NAME, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

N_max_test %>%
  ggplot() +
  geom_point(aes(x = species, y = N_prop, color = NA_L1NAME))+
  facet_wrap(~rank, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
  

