# Rarefaction analysis

library(dplyr)
library(tidyverse)
library(vegan)
library(mobr)
library(ggplot2)
library(betapart)
library(lme4)
library(performance)
library(tibble)
library(scales)

######### OBJECTIVE 1 - CONTINENTAL SCALE ########
######### ALL GIWs rarefaction ##############

wet_all <- readRDS("Intermediate_data/all_wet_check_for_analysis.RDS")

env <- readRDS("Data/earth_engine_env_data/env_matrix.RDS")

# need spatial coords for the LOCALITY_IDs
wet_coords <- readRDS("Intermediate_data/locality_ids_long_lat.RDS")

env <- env %>%
  left_join(wet_coords, by = "LOCALITY_ID")

## Creating a dataframe with just the LOCALITY_ID and ECOREGION
## to append to the wet_all df

ecoregion_loc <- env %>%
  select(LOCALITY_ID, NA_L1NAME)

wet_all <- wet_all %>%
  inner_join(ecoregion_loc, by = "LOCALITY_ID")


## Now a function that create comm 
## that can be used to plot the rarefaction curves
## but still matches the bootstrap we have by year

build_comm_mobr <- function(df, effort = 5, n_boot = 99, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  df %>%
    group_by(YEAR) %>%
    group_split() %>%
    map_dfr(function(year_dat) {
      year_val <- year_dat$YEAR[1]
      
      # Sites with enough checklists
      sites_ok <- year_dat %>%
        distinct(LOCALITY_ID, SAMPLING_EVENT_IDENTIFIER) %>%
        count(LOCALITY_ID, name = "n_checklists") %>%
        filter(n_checklists >= effort) %>%
        pull(LOCALITY_ID)
      
      year_dat <- year_dat %>% filter(LOCALITY_ID %in% sites_ok)
      if (n_distinct(year_dat$NA_L1NAME) < 2) return(NULL)
      
      map_dfr(seq_len(n_boot), function(b) {
        chosen_sites <- year_dat %>%
          distinct(NA_L1NAME, LOCALITY_ID) %>%
          group_by(NA_L1NAME) %>%
          slice_sample(n = 1) %>%
          pull(LOCALITY_ID)
        
        chosen_chk <- year_dat %>%
          filter(LOCALITY_ID %in% chosen_sites) %>%
          distinct(LOCALITY_ID, SAMPLING_EVENT_IDENTIFIER) %>%
          group_by(LOCALITY_ID) %>%
          slice_sample(n = effort) %>%
          ungroup()
        
        sampled <- year_dat %>%
          semi_join(chosen_chk, by = c("LOCALITY_ID", "SAMPLING_EVENT_IDENTIFIER"))
        
        comm <- sampled %>%
          group_by(LOCALITY_ID, COMMON_NAME) %>%
          summarise(abundance = sum(OBSERVATION_COUNT, na.rm = TRUE), .groups = "drop") %>%
          pivot_wider(names_from = COMMON_NAME, values_from = abundance, values_fill = 0)
        
        comm$YEAR <- year_val
        comm$bootstrap <- b
        comm
      })
    })
}

comm_df <- build_comm_mobr(wet_all)

# Still uneven by Year, since some years have more samples than others, let's downsample

comm_df_samp <- comm_df %>%
  group_by(YEAR) %>%
  slice_sample(n = min(table(comm_df$YEAR))) %>%
  ungroup()

# Now repalce all the NAs
comm_df_down <- comm_df_samp %>% 
  mutate(across(everything(), ~replace_na(.x, 0))) %>%
  mutate(sample_id = paste(LOCALITY_ID, YEAR, bootstrap, sep = "-"))

## What we want to do now, is just sample 
## 1 LOCALITY_ID out of this table

comm_df_uni <- comm_df_down %>%
  group_by(LOCALITY_ID) %>%
  slice_sample(n = 1) %>%
  ungroup

# make the rownames a unique sampling identifier
comm_df_final <- comm_df_uni %>%
  column_to_rownames(var = "sample_id") %>%
  select(-LOCALITY_ID, -YEAR, -bootstrap)

# now make the env df the same length
env_df <- comm_df_uni %>%
  select(LOCALITY_ID, YEAR, bootstrap, sample_id) %>%
  left_join(env, by = "LOCALITY_ID") %>%
  column_to_rownames(var = "sample_id") %>%
  mutate(dummy = "wetland")



mob_wet <- make_mob_in(comm_df_final, env_df, coord_names = c('LONGITUDE', 'LATITUDE'))

# Regular (non-spatial) SBR
rf_sbr <- rarefaction(mob_wet, method = "SBR")

# Spatial sample-based (kNN)
rf_ssbr_knn <- rarefaction(mob_wet, method = "sSBR", spat_algo = "kNN")

# Spatial sample-based (kNCN)
rf_ssbr_kncn <- rarefaction(mob_wet, method = "sSBR", spat_algo = "kNCN")

# Turn into data frames
df_sbr <- tibble(n = seq_along(rf_sbr), S = rf_sbr, method = "SBR")
df_ssbr_knn  <- tibble(n = seq_along(rf_ssbr_knn), S = rf_ssbr_knn, method = "sSBR_kNN")
df_ssbr_kncn <- tibble(n = seq_along(rf_ssbr_kncn), S = rf_ssbr_kncn, method = "sSBR_kNCN")

# Combine
df_all <- bind_rows(df_sbr, df_ssbr_knn, df_ssbr_kncn)

# Plot
ggplot(df_all, aes(x = n, y = S, color = method)) +
  geom_line(size = 1.1) +
  theme_minimal() +
  labs(
    x = "Number of samples",
    y = "Rarefied species richness",
    color = "Method",
    title = "Comparison: SBR, sSBR (kNN & kNCN)"
  ) +
  scale_color_manual(values = c(
    SBR = "black",
    sSBR_kNN = "blue",
    sSBR_kNCN = "red"
  )) +
  theme(legend.position = "bottom")


### now let's do this WITHIN ecoregion
## First let's filter down to our ecoregions that acutally worked for our analysis

env_df_count <- env_df %>%
  count(NA_L1NAME) %>%
  filter(n > 5)

env_df_eco <- env_df %>%
  filter(NA_L1NAME %in% env_df_count$NA_L1NAME)

comm_df_uni_eco <- comm_df_uni %>%
  filter(LOCALITY_ID %in% env_df_eco$LOCALITY_ID)


plot_rarefaction_ecoregion <- function(comm_df, env_df, ecoregion_name, coord_names = c("LONGITUDE", "LATITUDE")) {
  
  # Filter environment data for this ecoregion
  eco_env <- env_df %>%
    filter(NA_L1NAME == ecoregion_name)
  
  # Filter community matrix for the sites in this ecoregion
  comm_sub <- comm_df %>%
    filter(LOCALITY_ID %in% eco_env$LOCALITY_ID) %>%
    column_to_rownames(var = "sample_id") %>%
    select(-LOCALITY_ID, -YEAR, -bootstrap)
  
  # Build mobr object
  mob_in <- make_mob_in(comm_sub, eco_env, coord_names = coord_names)
  
  # Compute rarefaction curves
  rf_sbr      <- rarefaction(mob_in, method = "SBR")
  rf_ssbr_knn <- rarefaction(mob_in, method = "sSBR", spat_algo = "kNN")
  rf_ssbr_kncn <- rarefaction(mob_in, method = "sSBR", spat_algo = "kNCN")
  
  # Convert to tidy data frames
  df_sbr      <- tibble(n = seq_along(rf_sbr), S = rf_sbr, method = "SBR")
  df_ssbr_knn  <- tibble(n = seq_along(rf_ssbr_knn), S = rf_ssbr_knn, method = "sSBR_kNN")
  df_ssbr_kncn <- tibble(n = seq_along(rf_ssbr_kncn), S = rf_ssbr_kncn, method = "sSBR_kNCN")
  
  df_all <- bind_rows(df_sbr, df_ssbr_knn, df_ssbr_kncn)
  
  # Plot
  p <- ggplot(df_all, aes(x = n, y = S, color = method)) +
    geom_line(size = 1.1) +
    theme_minimal() +
    labs(
      x = "Number of samples",
      y = "Rarefied species richness",
      color = "Method",
      title = paste(ecoregion_name, ": SBR, sSBR (kNN & kNCN)")
    ) +
    scale_color_manual(values = c(
      SBR = "black",
      sSBR_kNN = "blue",
      sSBR_kNCN = "red"
    )) +
    theme(legend.position = "bottom")
  
  return(p)
}



unique_ecos <- unique(env_df_eco$NA_L1NAME)

plots <- map(unique_ecos, ~ plot_rarefaction_ecoregion(comm_df_uni_eco, env_df_eco, .x))

plots


## Additional function that truncates down to the minimum number of sites across ecoregions

plot_rarefaction_ecoregions <- function(comm_df, env_df, coord_names = c("LONGITUDE", "LATITUDE"), min_sites = 2) {
  
  # Compute minimum number of sites across ecoregions
  eco_nsites <- env_df %>%
    group_by(NA_L1NAME) %>%
    summarise(n_sites = n_distinct(LOCALITY_ID), .groups = "drop")
  
  # Minimum number of sites to truncate rarefaction curves
  min_n <- min(eco_nsites$n_sites)
  
  # Filter out ecoregions with too few sites
  valid_ecos <- eco_nsites %>% filter(n_sites >= min_sites) %>% pull(NA_L1NAME)
  
  # Loop over valid ecoregions
  rare_df <- map_dfr(valid_ecos, function(eco_name) {
    
    # Filter for this ecoregion
    eco_env <- env_df %>% filter(NA_L1NAME == eco_name)
    comm_sub <- comm_df %>%
      filter(LOCALITY_ID %in% eco_env$LOCALITY_ID) %>%
      column_to_rownames(var = "sample_id") %>%
      select(-LOCALITY_ID, -YEAR, -bootstrap)
    
    # Skip if still too few sites
    if(nrow(comm_sub) < min_sites) return(NULL)
    
    # Build mobr object
    mob_in <- make_mob_in(comm_sub, eco_env, coord_names = coord_names)
    
    # Safely compute rarefaction curves
    safe_raref <- function(method, spat_algo = NULL) {
      tryCatch(
        rarefaction(mob_in, method = method, spat_algo = spat_algo),
        error = function(e) NULL
      )
    }
    
    rf_sbr      <- safe_raref("SBR")
    rf_ssbr_knn <- safe_raref("sSBR", "kNN")
    rf_ssbr_kncn <- safe_raref("sSBR", "kNCN")
    
    # Truncate all curves to min_n
    rf_sbr      <- rf_sbr[1:min_n]
    rf_ssbr_knn <- rf_ssbr_knn[1:min_n]
    rf_ssbr_kncn <- rf_ssbr_kncn[1:min_n]
    
    # Combine into tidy df
    bind_rows(
      tibble(n = seq_along(rf_sbr), S = rf_sbr, method = "SBR", ecoregion = eco_name),
      tibble(n = seq_along(rf_ssbr_knn), S = rf_ssbr_knn, method = "sSBR_kNN", ecoregion = eco_name),
      tibble(n = seq_along(rf_ssbr_kncn), S = rf_ssbr_kncn, method = "sSBR_kNCN", ecoregion = eco_name)
    )
  })
  
  # Plot all ecoregions with facets
  ggplot(rare_df, aes(x = n, y = S, color = method)) +
    geom_line(size = 1) +
    facet_wrap(~ ecoregion, scales = "free_y") +
    theme_minimal(base_size = 14) +
    labs(
      x = "Number of samples",
      y = "Rarefied species richness",
      color = "Method",
      title = paste("Rarefaction Curves by Ecoregion (truncated to", min_n, "samples)")
    ) +
    scale_color_manual(values = c(
      SBR = "black",
      sSBR_kNN = "blue",
      sSBR_kNCN = "red"
    )) +
    theme(legend.position = "bottom")
}

plot_rarefaction_ecoregions(comm_df_uni_eco, env_df_eco, coord_names = c("LONGITUDE", "LATITUDE"))
