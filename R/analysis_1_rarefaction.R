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
library(patchwork)
library(sf)

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

# Now replace all the NAs
comm_df <- comm_df %>% 
  mutate(across(everything(), ~replace_na(.x, 0))) %>%
  mutate(sample_id = paste(LOCALITY_ID, YEAR, bootstrap, sep = "-"))

## What we want to do now, is just sample 
## 1 LOCALITY_ID out of this table
comm_df_uni <- comm_df %>%
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

env_df_filt <-env_df %>%
  select(LOCALITY_ID, LONGITUDE, LATITUDE, NA_L1NAME)


## Need to get the mean abundance from the community matrix
comm_sum <- rowSums(comm_df_final) %>%
  mean()


mob_wet <- make_mob_in(comm_df_final, env_df_filt, coord_names = c('LONGITUDE', 'LATITUDE'), 
                       latlong = TRUE)

# Spatial sample-based (kNN)
rf_ssbr_knn <- rarefaction(mob_wet, method = "sSBR", spat_algo = "kNN", latlong = TRUE)

# Random individual based rarefaction (IBR)
rf_ibr <- rarefaction(mob_wet, method = "IBR")

n_ibr <- 0:(length(rf_ibr)-1)

df_ibr <- tibble(n = n_ibr / comm_sum, S = rf_ibr)

# Turn into data frames
df_ssbr_knn  <- tibble(n = seq_along(rf_ssbr_knn), S = rf_ssbr_knn, method = "sSBR_kNN")


## Dan's suggested plot
#compare_samp_rarefaction(mob_wet)
#lines(1:length(rf_ibr) / comm_sum, rf_ibr, col='purple')

# Plot
ggplot(df_ssbr_knn, aes(x = n, y = S, color = method)) +
  geom_line(size = 1.1) +
  geom_line(data = df_ibr, aes(x = n, y = S, color = "IBR"), linewidth = 1.1) +
  scale_color_manual(values = c(
    "sSBR_kNN" = "black",
    "IBR" = "purple"
  )) +
  labs(
    x = "# of GIWs",
    y = "Species richness"
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "bottom")






# ---- BOOTSTRAP FUNCTION ----
run_single_boot <- function(comm_df, env, comm_sum) {
  
  # Select 1 row per LOCALITY_ID
  comm_df_uni <- comm_df %>%
    group_by(LOCALITY_ID) %>%
    slice_sample(n = 1) %>%
    ungroup()
  
  # Build comm matrix
  comm_df_final <- comm_df_uni %>%
    column_to_rownames(var = "sample_id") %>%
    select(-LOCALITY_ID, -YEAR, -bootstrap)
  
  ## Need to get the mean abundance from the community matrix
  comm_sum <- rowSums(comm_df_final) %>%
    mean()
  
  # Build env
  env_df <- comm_df_uni %>%
    select(LOCALITY_ID, YEAR, bootstrap, sample_id) %>%
    left_join(env, by = "LOCALITY_ID") %>%
    column_to_rownames(var = "sample_id") %>%
    mutate(dummy = "wetland")
  
  env_df_filt <- env_df %>%
    select(LOCALITY_ID, LONGITUDE, LATITUDE, NA_L1NAME)
  
  # Build mobr object
  mob_wet <- make_mob_in(
    comm_df_final, env_df_filt,
    coord_names = c("LONGITUDE", "LATITUDE"),
    latlong = TRUE
  )
  
  # ---- Rarefaction curves ----
  
  # Spatial SBR (kNN)
  rf_ssbr_knn <- rarefaction(mob_wet, method = "sSBR",
                             spat_algo = "kNN", latlong = TRUE)
  
  # IBR
  rf_ibr <- rarefaction(mob_wet, method = "IBR")
  n_ibr <- 1:(length(rf_ibr)-1)
  
  # Output tidy list
  list(
    ssbr = tibble(n = seq_along(rf_ssbr_knn), S = rf_ssbr_knn),
    ibr  = tibble(n = n_ibr / comm_sum, S = rf_ibr, pull = seq_along(rf_ibr))
  )
}

set.seed(123)

nboot <- 5

boot_results <- map(1:nboot, ~ run_single_boot(comm_df, env, comm_sum))

ssbr_all <- map_df(boot_results, "ssbr", .id = "rep")

ssbr_ci <- ssbr_all %>%
  group_by(n) %>%
  summarise(
    mean = mean(S),
    lower = quantile(S, 0.025),
    upper = quantile(S, 0.975)
  ) %>%
  mutate(method = "sSBR_kNN")

ibr_all <- map_df(boot_results, "ibr", .id = "rep")

# Find the minimum curve length across bootstrap runs
min_len <- ibr_all %>%
  group_by(rep) %>%
  summarise(L = max(pull)) %>%
  pull(L) %>%
  min()

# Truncate all IBR curves to that length
ibr_trunc <- ibr_all %>% 
  filter(pull <= min_len)

ibr_ci <- ibr_trunc %>%
  group_by(pull) %>%
  summarise(
    n = mean(n),
    mean = mean(S),
    lower = quantile(S, 0.025),
    upper = quantile(S, 0.975)
  ) %>%
  mutate(method = "IBR")




ggplot() +
  # sSBR CI band
  geom_ribbon(data = ssbr_ci,
              aes(x = n, ymin = lower, ymax = upper),
              fill = "grey80", alpha = 0.5) +
  geom_line(data = ssbr_ci, aes(x = n, y = mean), color = "black", linewidth = 1.1) +
  
  # IBR CI band
  geom_ribbon(data = ibr_ci,
              aes(x = n, ymin = lower, ymax = upper),
              fill = "purple", alpha = 0.2) +
  geom_line(data = ibr_ci, aes(x = n, y = mean), color = "purple", linewidth = 1.1) +
  
  labs(
    x = "# of GIWs",
    y = "Species richness"
  ) +
  theme_minimal(base_size = 10) +
  theme(legend.position = "none")




#### OBJECTIVE 2 - ECOREGION ########
### now let's do this WITHIN ecoregion
## First let's filter down to our ecoregions that actually worked for our analysis
env_df_count <- env_df_filt %>%
  count(NA_L1NAME) %>%
  filter(n > 5)

env_df_eco <- env_df_filt %>%
  filter(NA_L1NAME %in% env_df_count$NA_L1NAME)

comm_df_uni_eco <- comm_df_uni %>%
  filter(LOCALITY_ID %in% env_df_eco$LOCALITY_ID)

env_df_eco <- env_df_eco %>%
  mutate(
    LONGITUDE = as.numeric(LONGITUDE),
    LATITUDE  = as.numeric(LATITUDE)
  ) %>%
  filter(
    !is.na(LONGITUDE),
    !is.na(LATITUDE),
    LONGITUDE >= -180, LONGITUDE <= 180,
    LATITUDE >= -90, LATITUDE <= 90
  )


plot_rarefaction_ecoregion <- function(comm_df, env_df, ecoregion_name, coord_names = c("LONGITUDE", "LATITUDE")) {
  
  # Filter environment data for this ecoregion
  eco_env <- env_df %>%
    filter(NA_L1NAME == ecoregion_name)
  
  # Filter community matrix for the sites in this ecoregion
  comm_sub <- comm_df %>%
    filter(LOCALITY_ID %in% eco_env$LOCALITY_ID) %>%
    column_to_rownames(var = "sample_id") %>%
    select(-LOCALITY_ID, -YEAR, -bootstrap)
  
  #mean abundance per wetland site
  comm_mean <- rowSums(comm_sub) %>%
    mean()
  
  # Build mobr object
  mob_in <- make_mob_in(comm_sub, eco_env, coord_names = coord_names, latlong = TRUE)
  
  # Compute rarefaction curves
  rf_ssbr_knn <- rarefaction(mob_in, method = "sSBR", spat_algo = "kNN")

  # Prepend 0 so lines start at origin
  rf_ssbr_knn  <- c(0, rf_ssbr_knn)

  
  # Convert to tidy data frames
  df_ssbr_knn  <- tibble(n = seq_along(rf_ssbr_knn), S = rf_ssbr_knn, method = "sSBR_kNN")
  
  # curves for the IBR
  rf_ibr <- rarefaction(mob_in, method = "IBR")
  n_ibr <- 0:(length(rf_ibr)-1)
  df_ibr <- tibble(n = n_ibr / comm_mean, S = rf_ibr)
  
  
  ## Dan's suggested plot
  p <- ggplot(df_ssbr_knn, aes(x = n, y = S, color = method)) +
    geom_line(size = 1.1) +
    geom_line(data = df_ibr, aes(x = n, y = S, color = "IBR"), size = 1.1) +
    scale_color_manual(values = c(
      "sSBR_kNN" = "black",
      "IBR" = "purple"
    )) +
    labs(
      title = ecoregion_name,
      x = "# of GIWs",
      y = "Species richness"
    ) +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom",
          plot.title = element_text(size = 10))
  
  return(p)
}



unique_ecos <- unique(env_df_eco$NA_L1NAME)

plots <- map(unique_ecos, ~ plot_rarefaction_ecoregion(comm_df_uni_eco, env_df_eco, .x))

plots

(wrap_plots(plots[c(1,5,2, 6:7)]) +
    plot_layout(guides = "collect")) &
  theme(legend.position = "bottom")












## Additional function that truncates down to the minimum number of sites across ecoregions

plot_rarefaction_ecoregions <- function(comm_df, env_df, coord_names = c("LONGITUDE", "LATITUDE"), min_sites = 5) {
  # Compute number of sites per ecoregion
  eco_nsites <- env_df %>%
    group_by(NA_L1NAME) %>%
    summarise(n_sites = n_distinct(LOCALITY_ID), .groups = "drop")
  
  # Filter to ecoregions with at least 'min_sites'
  valid_ecos <- eco_nsites %>%
    filter(n_sites >= min_sites) %>%
    pull(NA_L1NAME)
  
  # Loop over ecoregions
  rare_plots <- map(valid_ecos, function(eco_name) {
    message("Processing ecoregion: ", eco_name)
    
    # Filter environmental and community data
    eco_env <- env_df %>% filter(NA_L1NAME == eco_name)
    comm_sub <- comm_df %>% filter(LOCALITY_ID %in% eco_env$LOCALITY_ID)
    
    comm_sub <- comm_df %>%
      filter(LOCALITY_ID %in% eco_env$LOCALITY_ID) %>%
      select(LOCALITY_ID, where(is.numeric))   # << NEW
    
    comm_mean <- comm_sub %>%
      select(where(is.numeric)) %>%           # << FIX
      rowSums() %>%
      mean()
    
    # Build mobr object
    mob_in <- make_mob_in(comm_sub, eco_env, coord_names = coord_names)
    
    # Compute rarefaction curves
    rf_sbr       <- rarefaction(mob_in, method = "SBR")
    rf_ssbr_knn  <- rarefaction(mob_in, method = "sSBR", spat_algo = "kNN")
    rf_ssbr_kncn <- rarefaction(mob_in, method = "sSBR", spat_algo = "kNCN")
    rf_ibr       <- rarefaction(mob_in, method = "IBR")
    
    # Truncate all curves to match shortest
    min_len <- min(length(rf_sbr), length(rf_ssbr_knn), length(rf_ssbr_kncn), length(rf_ibr))
    rf_sbr       <- rf_sbr[1:min_len]
    rf_ssbr_knn  <- rf_ssbr_knn[1:min_len]
    rf_ssbr_kncn <- rf_ssbr_kncn[1:min_len]
    rf_ibr       <- rf_ibr[1:min_len]
    
    # Prepend 0 so lines start at origin
    rf_sbr       <- c(0, rf_sbr)
    rf_ssbr_knn  <- c(0, rf_ssbr_knn)
    rf_ssbr_kncn <- c(0, rf_ssbr_kncn)

    
    # Build x-values (sampling effort)
    n <- 0:(min_len)  # same for SBR, kNN, kNCN
    n_ibr <- n / comm_mean  # x-axis scaled for IBR
    
    # Build tidy dataframe
    df_all <- bind_rows(
      tibble(n = n,     S = rf_sbr,       method = "SBR"),
      tibble(n = n,     S = rf_ssbr_knn,  method = "sSBR_kNN"),
      tibble(n = n,     S = rf_ssbr_kncn, method = "sSBR_kNCN"),
      tibble(n = n_ibr, S = rf_ibr,       method = "IBR")
    )
    
    # Plot
    ggplot(df_all, aes(x = n, y = S, color = method)) +
      geom_line(size = 1.1) +
      scale_color_manual(values = c(
        "SBR" = "black",
        "sSBR_kNN" = "red",
        "sSBR_kNCN" = "blue",
        "IBR" = "purple"
      )) +
      labs(
        title = eco_name,
        x = "Sampling effort",
        y = "Species richness"
      ) +
      theme_minimal(base_size = 10) +
      theme(legend.position = "bottom")
  })
  
  return(rare_plots)
}




plots <- plot_rarefaction_ecoregions(comm_df_uni_eco, env_df_eco, coord_names = c("LONGITUDE", "LATITUDE"))

(wrap_plots(plots[c(1:3, 5, 7)]) +
  plot_layout(guides = "collect")) &
  theme(legend.position = "bottom")
