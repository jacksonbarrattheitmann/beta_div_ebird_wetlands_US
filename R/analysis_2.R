library(mobr)
library(dplyr)
library(tidyr)
library(geosphere)
library(purrr)
library(sf)
library(units)

####### OBJECTIVE 2 - ECOREGION scale ##########
wet_eco <- readRDS("Intermediate_data/wet_comm_ecoregion_10_summarized.RDS") %>%
  column_to_rownames(var = "LOCALITY_ID")

# the ENV data with ECOREGION = NA_L1NAME column
env_eco <- readRDS("Intermediate_data/env_ecoregions_10_GIWs.RDS") 

# need spatial coords for the LOCALITY_IDs
wet_coords <- readRDS("Intermediate_data/locality_ids_long_lat.RDS")

# should join the coords with the env dataframe
env_eco <- env_eco %>%
  inner_join(wet_coords, by = "LOCALITY_ID")

## Now the function

calculate_beta_ecoregion_sf <- function(wet_eco, env_eco, n_reps = 999, effort = 25, C_target_gamma = 0.75) {
  
  env_eco <- env_eco %>%
    filter(LOCALITY_ID %in% rownames(wet_eco)) %>%
    filter(!is.na(LATITUDE), !is.na(LONGITUDE)) %>%
    filter(LATITUDE >= -90, LATITUDE <= 90,
           LONGITUDE >= -180, LONGITUDE <= 180)
  
  ecoregions <- unique(env_eco$NA_L1NAME)
  results_df <- data.frame()
  
  for (eco in ecoregions) {
    message("Processing ecoregion: ", eco)
    
    env_sub <- env_eco %>% filter(NA_L1NAME == eco)
    if (nrow(env_sub) < 5) next
    
    env_sub_sf <- st_as_sf(env_sub, coords = c("LONGITUDE", "LATITUDE"), crs = 4326) %>%
      st_transform(crs = 5070)
    
    dist_matrix <- st_distance(env_sub_sf)
    mean_dist <- mean(dist_matrix[upper.tri(dist_matrix)], na.rm = TRUE)
    
    species_sub <- wet_eco[rownames(wet_eco) %in% env_sub$LOCALITY_ID, , drop = FALSE]
    
    replicate_counter <- 0
    max_anchor_attempts <- 1000
    
    while (replicate_counter < n_reps) {
      anchor_attempts <- 0
      valid <- FALSE
      
      while (!valid && anchor_attempts < max_anchor_attempts) {
        anchor_idx <- sample(1:nrow(env_sub_sf), 1)
        dists <- st_distance(env_sub_sf[anchor_idx, ], env_sub_sf)
        candidates <- which(dists <= mean_dist & dists > set_units(0, "m"))
        
        if (length(candidates) >= 4) {
          selected_idxs <- c(anchor_idx, sample(candidates, 4))
          valid <- TRUE
        }
        
        anchor_attempts <- anchor_attempts + 1
      }
      
      if (!valid) break
      
      selected_ids <- env_sub$LOCALITY_ID[selected_idxs]
      selected_comm <- wet_eco[rownames(wet_eco) %in% selected_ids, , drop = FALSE]
      selected_comm <- selected_comm[, colSums(selected_comm) > 0, drop = FALSE]
      
      if (ncol(selected_comm) == 0 || nrow(selected_comm) != 5) next
      
      selected_comm <- as.data.frame(selected_comm)
      selected_comm[] <- lapply(selected_comm, as.numeric)
      rownames(selected_comm) <- selected_ids
      
      if (any(is.na(selected_comm)) || any(rowSums(selected_comm) == 0)) next
      
      comm_div <- tryCatch({
        calc_comm_div(
          selected_comm,
          index = c('S', 'S_n', 'S_PIE', 'S_C'),
          extrapolate = TRUE,
          effort = 25,
          scales = c("beta"),
          C_target_gamma = C_target_gamma
        )
      }, error = function(e) NULL)
      
      if (is.null(comm_div) || nrow(comm_div) == 0) next
      
      comm_div$replicate <- replicate_counter + 1
      comm_div$NA_L1NAME <- eco
      
      results_df <- bind_rows(results_df, comm_div)
      replicate_counter <- replicate_counter + 1
    }
    
    if (replicate_counter < n_reps) {
      warning(paste("Only", replicate_counter, "valid replicates for", eco))
    }
  }
  
  return(results_df)
}

beta_results <- calculate_beta_ecoregion_sf(wet_eco, env_eco, n_reps = 100)

# error bars confidence intervals 95%
beta_div_error <- beta_results %>%
  group_by(index, NA_L1NAME) %>%
  summarize(
    mean = mean(value, na.rm = TRUE),
    lower = quantile(value, 0.025, na.rm = TRUE),
    upper = quantile(value, 0.975, na.rm = TRUE)
  )

ggplot() +
  geom_jitter(data = beta_results, aes(x = index, y = value, color = NA_L1NAME)) +
  geom_point(data = beta_div_error, aes(x = index, y = mean), color = "black")+
  geom_errorbar(data = beta_div_error, aes(x = index, ymin = lower, ymax = upper, width = 0.2), color = "black", alpha = 0.75, linewidth = 1) +
  geom_hline(yintercept = 1, color = "darkred", linetype = "dashed", linewidth = 1) +
  facet_wrap(~NA_L1NAME) +
  theme_bw() +
  theme(legend.position = "none")  +
  ylab("Value") +
  xlab("Beta Diversity Index") +
  scale_x_discrete(labels = c(
    "beta_S" = "βS",
    "beta_S_n" = "βSn",
    "beta_S_PIE" = "βSPIE", 
    "beta_S_C" = "βC"))

ggsave("FigX_betas_ECOREGION_scale.png", width = 8, height = 6,
       bg = "transparent")

## the model

for (idx in unique(beta_results$index)) {
  cat("\n===== Index:", idx, "=====\n")
  
  df <- beta_results %>% filter(index == idx)
  
  aov_model <- aov(value ~ NA_L1NAME, data = df)
  print(summary(aov_model))
  
  tukey <- TukeyHSD(aov_model)
  print(tukey)
  
  # Plot with title
  plot(tukey, las = 1)
  title(main = paste("Tukey HSD –", idx))
}

   