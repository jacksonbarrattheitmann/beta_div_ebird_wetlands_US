# analysis script

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
######### ALL GIWs beta ##############

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
  

##### FUNCTION for calculating beta 

calc_beta_by_year <- function(df, effort = 5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)  # reproducibility if seed supplied
  
  df %>%
    group_by(YEAR) %>%
    group_split() %>%
    map_dfr(function(year_dat) {
      year_val <- year_dat$YEAR[1]
      
      # Sites with at least `effort` unique checklists
      sites_ok <- year_dat %>%
        distinct(LOCALITY_ID, SAMPLING_EVENT_IDENTIFIER) %>%
        count(LOCALITY_ID, name = "n_checklists") %>%
        filter(n_checklists >= effort) %>%
        pull(LOCALITY_ID)
      
      if (length(sites_ok) < 2) return(NULL)  # need ≥2 sites
      
      # Sample effort checklists per site
      chosen_chk <- year_dat %>%
        filter(LOCALITY_ID %in% sites_ok) %>%
        distinct(LOCALITY_ID, SAMPLING_EVENT_IDENTIFIER) %>%
        group_by(LOCALITY_ID) %>%
        slice_sample(n = effort) %>%
        ungroup()
      
      sampled <- year_dat %>%
        semi_join(chosen_chk, by = c("LOCALITY_ID", "SAMPLING_EVENT_IDENTIFIER"))
      
      # Collapse 5 checklists into a single row per site
      comm <- sampled %>%
        group_by(LOCALITY_ID, COMMON_NAME) %>%
        summarise(abundance = sum(OBSERVATION_COUNT, na.rm = TRUE), .groups = "drop") %>%
        pivot_wider(names_from = COMMON_NAME, values_from = abundance, values_fill = 0) %>%
        as.data.frame()
      
      rownames(comm) <- comm$LOCALITY_ID
      comm$LOCALITY_ID <- NULL  # now rows = sites, cols = species
      
      # Beta diversity
      beta_out <- calc_comm_div(
        comm,
        index = c("S", "S_n", "S_PIE", "S_C"),
        extrapolate = TRUE,
        effort = 25,
        scales = "beta",
        C_target_gamma = 0.75
      )
      
      tibble(
        YEAR = year_val,
        n_sites = length(sites_ok),
        scale = beta_out$scale,
        index = beta_out$index,
        sample_size = beta_out$sample_size,
        effort = beta_out$effort,
        gamma_coverage = beta_out$gamma_coverage,
        value = beta_out$value
      )
    })
}


betas_ALL <- calc_beta_by_year(wet_all, effort = 5, seed = 5)

########## PLOTS for OBJECTIVE 1 ##################


ggplot(data = betas_ALL) +
  geom_point(aes(x = index, y = value)) + theme_bw() +
  geom_hline(yintercept = 1, color = "darkred", linetype = "dashed", linewidth = 1) +
  xlab("Beta Diversity Index") +
  ylab("Value") +
  scale_x_discrete(  labels = c(
    "beta_S" = "βS",
    "beta_S_n" = "βSn",
    "beta_S_PIE" = "βSPIE", 
    "beta_S_C" = "βC"
  )) +
  theme(legend.position = "none",
        axis.title.x = element_text(
          margin = margin(t = 15)),
        axis.title.y = element_text(
          margin = margin(r = 15)))


ggsave("Fig1_betas_CONTINENTAL_SCALE.png", width = 6, height = 4,
       bg = "transparent")


### Let's do one more sensitivity analysis
### to deal with the distance-decay function, so let's just compare 1 wetland from each Ecoregion
### and bootstrap this a bunch of times


calc_beta_by_year_boot <- function(df, effort = 5, n_boot = 99, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)  # reproducibility
  
  df %>%
    group_by(YEAR) %>%
    group_split() %>%
    map_dfr(function(year_dat) {
      year_val <- year_dat$YEAR[1]
      
      # Sites with at least `effort` checklists
      sites_ok <- year_dat %>%
        distinct(LOCALITY_ID, SAMPLING_EVENT_IDENTIFIER) %>%
        count(LOCALITY_ID, name = "n_checklists") %>%
        filter(n_checklists >= effort) %>%
        pull(LOCALITY_ID)
      
      # Keep only eligible sites
      year_dat <- year_dat %>%
        filter(LOCALITY_ID %in% sites_ok)
      
      if (n_distinct(year_dat$NA_L1NAME) < 2) return(NULL) # need ≥2 ecoregions
      
      # Bootstrap loop
      map_dfr(seq_len(n_boot), function(b) {
        
        # Select one LOCALITY_ID per ecoregion
        chosen_sites <- year_dat %>%
          distinct(NA_L1NAME, LOCALITY_ID) %>%
          group_by(NA_L1NAME) %>%
          slice_sample(n = 1) %>%
          pull(LOCALITY_ID)
        
        # Sample `effort` checklists from each chosen site
        chosen_chk <- year_dat %>%
          filter(LOCALITY_ID %in% chosen_sites) %>%
          distinct(LOCALITY_ID, SAMPLING_EVENT_IDENTIFIER) %>%
          group_by(LOCALITY_ID) %>%
          slice_sample(n = effort) %>%
          ungroup()
        
        sampled <- year_dat %>%
          semi_join(chosen_chk, by = c("LOCALITY_ID", "SAMPLING_EVENT_IDENTIFIER"))
        
        # Collapse 5 checklists into one row per site
        comm <- sampled %>%
          group_by(LOCALITY_ID, COMMON_NAME) %>%
          summarise(abundance = sum(OBSERVATION_COUNT, na.rm = TRUE), .groups = "drop") %>%
          pivot_wider(names_from = COMMON_NAME, values_from = abundance, values_fill = 0) %>%
          as.data.frame()
        
        rownames(comm) <- comm$LOCALITY_ID
        comm$LOCALITY_ID <- NULL
        
        # Beta diversity
        beta_out <- calc_comm_div(
          comm,
          index = c("S", "S_n", "S_PIE", "S_C"),
          extrapolate = TRUE,
          effort = 25,
          scales = "beta",
          C_target_gamma = 0.75
        )
        
        tibble(
          YEAR = year_val,
          bootstrap = b,
          n_sites = length(chosen_sites),
          scale = beta_out$scale,
          index = beta_out$index,
          sample_size = beta_out$sample_size,
          effort = beta_out$effort,
          gamma_coverage = beta_out$gamma_coverage,
          value = beta_out$value
        )
      })
    })
}

betas_boot <- calc_beta_by_year_boot(wet_all, effort = 5, n_boot = 25, seed = NULL)



# calculating the means and error bars for plotting
wet_div_error <- betas_boot %>%
  group_by(index) %>%
  summarize(
    mean = mean(value, na.rm = TRUE),
    lower = quantile(value, 0.025, na.rm = TRUE),
    upper = quantile(value, 0.975, na.rm = TRUE)
  )


# Plotting the results
ggplot() +
  geom_jitter(data = betas_boot, aes(x = index, y = value), alpha = 0.25, width = 0.2) +
  geom_point(data = wet_div_error, aes(x = index, y = mean), color = "darkred", size = 3, alpha = 1) +
  geom_errorbar(data = wet_div_error, aes(x = index, ymin = lower, ymax = upper, width = 0.2), color = "darkred", linewidth = 1) +
  geom_hline(yintercept = 1, color = "dodgerblue", linetype = "dashed", linewidth = 1) +
  theme_bw() +
  theme(legend.position = "none") +
  ylab("Value") +
  xlab("Diversity Index") +
  scale_x_discrete(labels = c(
    "beta_S" = "βS",
    "beta_S_n" = "βSn",
    "beta_S_PIE" = "βSPIE", 
    "beta_S_C" = "βC"))



######## SUPPLEMENTARY FIGURES ############


##### Let's gut check by plotting by YEAR for each Index, to see if anything's changing

ggplot(data = betas_ALL) +
  geom_point(aes(x = YEAR, y = value)) +
  facet_wrap(~index) +
  theme_bw() +
  geom_hline(yintercept = 1, color = "darkred", linetype = "dashed", linewidth = 1) +
  xlab("Year") +
  ylab("Value")

# really looks the same across years, even as we accumulate more samples from other ecoregions
# i.e., we are not seeing high turnover between wetlands due to aggregation changes
# mostly due to changes in the abundance of common species



##### Let's gut check by plotting by YEAR for each Index, to see if anything's changing

ggplot(data = betas_boot) +
  geom_point(aes(x = YEAR, y = value)) +
  facet_wrap(~index) +
  theme_bw() +
  geom_hline(yintercept = 1, color = "darkred", linetype = "dashed", linewidth = 1) +
  xlab("Year") +
  ylab("Value")


#### Figure of sample sites per YEAR

ggplot(data = betas_ALL) +
  geom_point(aes(x = YEAR, y = n_sites)) +
  theme_bw() +
  xlab("Year") +
  ylab("Number of GIWs used to calculate β")


#### SENSITIVITY TEST for Sn ######

calc_beta_by_year_sn <- function(df, effort = 5, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)  # reproducibility if seed supplied
  
  df %>%
    group_by(YEAR) %>%
    group_split() %>%
    map_dfr(function(year_dat) {
      year_val <- year_dat$YEAR[1]
      
      # Sites with at least `effort` unique checklists
      sites_ok <- year_dat %>%
        distinct(LOCALITY_ID, SAMPLING_EVENT_IDENTIFIER) %>%
        count(LOCALITY_ID, name = "n_checklists") %>%
        filter(n_checklists >= effort) %>%
        pull(LOCALITY_ID)
      
      if (length(sites_ok) < 2) return(NULL)  # need ≥2 sites
      
      # Sample effort checklists per site
      chosen_chk <- year_dat %>%
        filter(LOCALITY_ID %in% sites_ok) %>%
        distinct(LOCALITY_ID, SAMPLING_EVENT_IDENTIFIER) %>%
        group_by(LOCALITY_ID) %>%
        slice_sample(n = effort) %>%
        ungroup()
      
      sampled <- year_dat %>%
        semi_join(chosen_chk, by = c("LOCALITY_ID", "SAMPLING_EVENT_IDENTIFIER"))
      
      # Collapse 5 checklists into a single row per site
      comm <- sampled %>%
        group_by(LOCALITY_ID, COMMON_NAME) %>%
        summarise(abundance = sum(OBSERVATION_COUNT, na.rm = TRUE), .groups = "drop") %>%
        pivot_wider(names_from = COMMON_NAME, values_from = abundance, values_fill = 0) %>%
        as.data.frame()
      
      rownames(comm) <- comm$LOCALITY_ID
      comm$LOCALITY_ID <- NULL  # now rows = sites, cols = species
      
      # Beta diversity
      beta_out <- calc_comm_div(
        comm,
        index = c("S_n"),
        extrapolate = TRUE,
        effort = c(10, 25, 50),
        scales = "beta",
        C_target_gamma = 0.75
      )
      
      tibble(
        YEAR = year_val,
        n_sites = length(sites_ok),
        scale = beta_out$scale,
        index = beta_out$index,
        sample_size = beta_out$sample_size,
        effort = beta_out$effort,
        gamma_coverage = beta_out$gamma_coverage,
        value = beta_out$value
      )
    })
}


betas_sn_testing <- calc_beta_by_year_sn(wet_all, seed = 5)

# calculating the means and error bars for plotting
wet_div_error <- betas_sn_testing %>%
  group_by(index, effort) %>%
  summarize(
    mean = mean(value, na.rm = TRUE),
    lower = quantile(value, 0.025, na.rm = TRUE),
    upper = quantile(value, 0.975, na.rm = TRUE)
  )


# Plotting the results
ggplot() +
  geom_jitter(data = betas_sn_testing, aes(x = index, y = value), alpha = 0.25, width = 0.2) +
  geom_point(data = wet_div_error, aes(x = index, y = mean), color = "darkred", size = 3, alpha = 1) +
  geom_errorbar(data = wet_div_error, aes(x = index, ymin = lower, ymax = upper, width = 0.2), color = "darkred", linewidth = 1) +
  geom_hline(yintercept = 1, color = "dodgerblue", linetype = "dashed", linewidth = 1) +
  facet_wrap(~effort) +
  theme_bw() +
  theme(legend.position = "none") +
  ylab("Value") +
  xlab("Diversity Index") +
  scale_x_discrete(labels = c(
    "beta_S_n" = "βSn"))



######## RAREFACTION ###############
