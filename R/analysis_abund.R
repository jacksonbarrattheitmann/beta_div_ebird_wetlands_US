library(dplyr)
library(tidyverse)
library(ggplot2)

### Checking out the Abundance data

wet_eco <- readRDS("Intermediate_data/wet_comm_ecoregion_10_summarized.RDS") %>%
  column_to_rownames(var = "LOCALITY_ID")

# the ENV data with ECOREGION = NA_L1NAME column
env_eco <- readRDS("Intermediate_data/env_ecoregions_10_GIWs.RDS") 


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
  mutate(rank_score = 6 - rank) %>%  # 6 - rank gives the scoring scheme you described
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
  

