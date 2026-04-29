library(mobr)
library(dplyr)
library(tidyr)
library(geosphere)
library(purrr)
library(sf)
library(units)
library(ggplot2)
library(tibble)
library(multcompView)
library(broom)

####### OBJECTIVE 2 - ECOREGION scale ##########
wet_eco_summer <- readRDS("Intermediate_data/all_wet_check_for_analysis_SUMMER.RDS")

wet_eco_fall <- readRDS("Intermediate_data/all_wet_check_for_analysis_FALL.RDS")

wet_eco_spring <- readRDS("Intermediate_data/all_wet_check_for_analysis_SPRING.RDS")

wet_eco_winter <- readRDS("Intermediate_data/all_wet_check_for_analysis_WINTER.RDS")

# the ENV data with ECOREGION = NA_L1NAME column
env <- readRDS("Data/old_data/earth_engine_env_data/env_matrix.RDS")

# need spatial coords for the LOCALITY_IDs
wet_coords <- readRDS("Intermediate_data/locality_ids_long_lat.RDS")

# should join the coords with the env dataframe
env_eco <- env %>%
  inner_join(wet_coords, by = "LOCALITY_ID")


wet_eco_summer <- wet_eco_summer %>%
  inner_join(env_eco, by = "LOCALITY_ID")

wet_eco_fall <- wet_eco_fall %>%
  inner_join(env_eco, by = "LOCALITY_ID")

wet_eco_spring <- wet_eco_spring %>%
  inner_join(env_eco, by = "LOCALITY_ID")

wet_eco_winter <- wet_eco_winter %>%
  inner_join(env_eco, by = "LOCALITY_ID") 

############# FUNCTION ##################

calculate_beta_ecoregion_sf <- function(wet_long,
                                        env_eco,
                                        n_reps = 999,
                                        C_target_gamma = 0.75,
                                        max_anchor_attempts = 1000) {
  
  results_df <- tibble()
  
  # Clean env metadata to only LOCALITY_IDs present in wet_long
  env_eco <- env_eco %>%
    filter(LOCALITY_ID %in% wet_long$LOCALITY_ID) %>%
    filter(!is.na(LATITUDE), !is.na(LONGITUDE)) %>%
    filter(LATITUDE >= -90, LATITUDE <= 90,
           LONGITUDE >= -180, LONGITUDE <= 180)
  
  ecoregions <- unique(env_eco$NA_L1NAME)
  
  for (eco in ecoregions) {
    message("Processing ecoregion: ", eco)
    env_sub <- env_eco %>% filter(NA_L1NAME == eco)
    if (nrow(env_sub) < 5) next
    
    env_sub_sf <- st_as_sf(env_sub, coords = c("LONGITUDE", "LATITUDE"), crs = 4326) %>%
      st_transform(crs = 5070)
    
    # Mean pairwise distance per ecoregion
    dist_matrix <- st_distance(env_sub_sf)
    mean_dist_num <- mean(as.numeric(dist_matrix[upper.tri(dist_matrix)]), na.rm = TRUE)
    
    replicate_counter <- 0
    
    while (replicate_counter < n_reps) {
      anchor_attempts <- 0
      valid <- FALSE
      
      while (!valid && anchor_attempts < max_anchor_attempts) {
        # 1) Pick an anchor LOCALITY_ID
        anchor_idx <- sample(seq_len(nrow(env_sub_sf)), 1)
        anchor_id <- env_sub$LOCALITY_ID[anchor_idx]
        
        # 2) Sample a YEAR for anchor LOCALITY_ID (pre-filtered, guaranteed ≥5 events)
        valid_years <- wet_long %>%
          filter(LOCALITY_ID == anchor_id) %>%
          pull(YEAR) %>%
          unique()
        
        if (length(valid_years) == 0) {
          anchor_attempts <- anchor_attempts + 1
          next
        }
        anchor_year <- sample(valid_years, 1)
        
        # 3) Find candidate nearby LOCALITY_IDs within mean distance
        dists <- st_distance(env_sub_sf[anchor_idx, ], env_sub_sf)
        dists_num <- as.numeric(dists)
        candidates <- which(dists_num <= mean_dist_num & dists_num > 0)
        
        # 4) Keep only candidates with ≥5 events in anchor_year
        # ensure candidates is a plain integer vector
        candidates <- unlist(candidates, use.names = FALSE)
        
        candidate_idxs_filtered <- candidates[
          vapply(candidates, function(i) {
            lid <- env_sub$LOCALITY_ID[i]
            n_events <- wet_long %>%
              filter(LOCALITY_ID == lid, YEAR == anchor_year) %>%
              distinct(SAMPLING_EVENT_IDENTIFIER) %>%
              nrow()
            n_events >= 5
          }, logical(1))
        ]
        
        if (length(candidate_idxs_filtered) < 4) {
          anchor_attempts <- anchor_attempts + 1
          next
        }
        
        # 5) Select 4 neighbors
        other_idxs <- sample(candidate_idxs_filtered, 4)
        selected_idxs <- c(anchor_idx, other_idxs)
        selected_ids <- env_sub$LOCALITY_ID[selected_idxs]
        
        # 6) Sample 5 events per LOCALITY_ID
        events_by_loc <- purrr::map(selected_ids, function(lid) {
          x <- wet_long %>%
            filter(LOCALITY_ID == lid, YEAR == anchor_year) %>%
            distinct(SAMPLING_EVENT_IDENTIFIER) %>%
            pull(SAMPLING_EVENT_IDENTIFIER)
          
          if (length(x) < 5) return(NULL)   # skip if not enough events
          sample(x, size = 5)
        })
        
        # Skip replicate if any site has <5 events
        if (any(sapply(events_by_loc, is.null))) {
          anchor_attempts <- anchor_attempts + 1
          next
        }
        
        events_map <- tibble(
          LOCALITY_ID = rep(selected_ids, each = 5),
          SAMPLING_EVENT_IDENTIFIER = unlist(events_by_loc)
        )
        
        # 7) Collapse each LOCALITY_ID into a single row
        comm_matrix_df <- wet_long %>%
          filter(YEAR == anchor_year,
                 SAMPLING_EVENT_IDENTIFIER %in% events_map$SAMPLING_EVENT_IDENTIFIER,
                 LOCALITY_ID %in% selected_ids) %>%
          group_by(LOCALITY_ID, COMMON_NAME) %>%
          summarise(abundance = sum(OBSERVATION_COUNT, na.rm = TRUE), .groups = "drop") %>%
          pivot_wider(names_from = COMMON_NAME,
                      values_from = abundance,
                      values_fill = list(abundance = 0))
        
        if (!"LOCALITY_ID" %in% names(comm_matrix_df)) {
          anchor_attempts <- anchor_attempts + 1
          next
        }
        
        comm_mat <- comm_matrix_df %>%
          column_to_rownames("LOCALITY_ID") %>%
          as.data.frame(stringsAsFactors = FALSE)
        
        # Drop all-zero species columns
        if (ncol(comm_mat) > 0) {
          comm_mat <- comm_mat[, colSums(comm_mat) > 0, drop = FALSE]
        }
        
        # Must be exactly 5 localities and non-empty rows
        if (nrow(comm_mat) == 5 &&
            ncol(comm_mat) > 0 &&
            all(rowSums(comm_mat) > 0)) {
          valid <- TRUE
        } else {
          anchor_attempts <- anchor_attempts + 1
        }
      } # end anchor_attempts
      
      if (!valid) break
      
      # 8) Compute diversity using MOBR
      comm_div <- tryCatch({
        calc_comm_div(
          comm_mat,
          index = c("S", "S_PIE", "S_C"),
          extrapolate = TRUE,
          scales = c("beta"),
          C_target_gamma = C_target_gamma
        )
      }, error = function(e) NULL)
      
      if (is.null(comm_div) || nrow(comm_div) == 0) {
        next
      }
      
      comm_div <- as_tibble(comm_div)
      comm_div$replicate <- replicate_counter + 1
      comm_div$NA_L1NAME <- eco
      comm_div$YEAR <- anchor_year
      comm_div$ANCHOR_ID <- anchor_id
      
      results_df <- bind_rows(results_df, comm_div)
      replicate_counter <- replicate_counter + 1
    } # end while replicates
    
    if (replicate_counter < n_reps) {
      warning(sprintf("Only %d valid replicates for ecoregion %s", replicate_counter, eco))
    }
  } # end ecoregions
  
  return(results_df)
}


# SUMMER
beta_results_summer <- calculate_beta_ecoregion_sf(wet_eco_summer, env_eco, n_reps = 99, max_anchor_attempts = 99)

beta_results_summer <- beta_results_summer %>%
  mutate(SEASON = "Summer")

# FALL
beta_results_fall <- calculate_beta_ecoregion_sf(wet_eco_fall, env_eco, n_reps = 99, max_anchor_attempts = 99)

beta_results_fall <- beta_results_fall %>%
  mutate(SEASON = "Fall")

# SPRING
beta_results_spring <- calculate_beta_ecoregion_sf(wet_eco_spring, env_eco, n_reps = 99, max_anchor_attempts = 99)

beta_results_spring <- beta_results_spring %>%
  mutate(SEASON = "Spring")

# WInter
beta_results_winter <- calculate_beta_ecoregion_sf(wet_eco_winter, env_eco, n_reps = 99, max_anchor_attempts = 99)

beta_results_winter <- beta_results_winter %>%
  mutate(SEASON = "Winter")

# Combine

beta_results_eco <- rbind(beta_results_winter, beta_results_spring, beta_results_summer, beta_results_fall)

#saveRDS(beta_results_eco, "Intermediate_data/beta_results_eco.RDS")

beta_results <- readRDS("Intermediate_data/beta_results_eco.RDS") %>%
  filter(SEASON == "Summer")

# error bars confidence intervals 95%
beta_div_error <- beta_results %>%
  group_by(index, NA_L1NAME, SEASON) %>%
  summarize(
    mean = mean(value, na.rm = TRUE),
    lower = quantile(value, 0.025, na.rm = TRUE),
    upper = quantile(value, 0.975, na.rm = TRUE)
  )


ggplot() +
  geom_jitter(
    data = beta_results,
    aes(
      x = factor(index, levels = c("beta_S", "beta_S_PIE", "beta_S_C")),
      y = value,
      color = NA_L1NAME
    ),
    alpha = 0.1,
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.5)
  ) +
  geom_point(
    data = beta_div_error,
    aes(
      x = factor(index, levels = c("beta_S", "beta_S_PIE", "beta_S_C")),
      y = mean,
      color = NA_L1NAME
    ),
    size = 3,
    shape = 21,
    fill = "white",
    stroke = 1.2,
    position = position_dodge(width = 0.5)
  ) +
  geom_errorbar(
    data = beta_div_error,
    aes(
      x = factor(index, levels = c("beta_S", "beta_S_PIE", "beta_S_C")),
      ymin = lower,
      ymax = upper,
      color = NA_L1NAME
    ),
    linewidth = 1.2,
    width = 0.2,
    position = position_dodge(width = 0.5)
  ) +
  geom_hline(yintercept = 1, color = "dodgerblue", linetype = "dashed", linewidth = 1) +
  facet_wrap(~NA_L1NAME) +
  theme_bw() +
  ylab("Value") +
  xlab("Diversity Index") +
  scale_x_discrete(labels = c(
    "beta_S" = "βS",
    "beta_S_PIE" = "βSPIE",
    "beta_S_C" = "βC"
  )) +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "none")







beta_results <- readRDS("Intermediate_data/beta_results_eco.RDS")

# remove MED CALI because we didn't have enough data in summer months

beta_results <- beta_results %>%
  filter(!NA_L1NAME == "MEDITERRANEAN CALIFORNIA")

# error bars confidence intervals 95%
beta_div_error <- beta_results %>%
  group_by(index, NA_L1NAME, SEASON) %>%
  summarize(
    mean = mean(value, na.rm = TRUE),
    lower = quantile(value, 0.025, na.rm = TRUE),
    upper = quantile(value, 0.975, na.rm = TRUE)
  )


ggplot() +
  geom_jitter(
    data = beta_results,
    aes(
      x = factor(index, levels = c("beta_S", "beta_S_PIE", "beta_S_C")),
      y = value,
      color = SEASON
    ),
    alpha = 0.1,
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.5)
  ) +
  geom_point(
    data = beta_div_error,
    aes(
      x = factor(index, levels = c("beta_S", "beta_S_PIE", "beta_S_C")),
      y = mean,
      color = SEASON
    ),
    size = 3,
    shape = 21,
    fill = "white",
    stroke = 1.2,
    position = position_dodge(width = 0.5)
  ) +
  geom_errorbar(
    data = beta_div_error,
    aes(
      x = factor(index, levels = c("beta_S", "beta_S_PIE", "beta_S_C")),
      ymin = lower,
      ymax = upper,
      color = SEASON
    ),
    linewidth = 1.2,
    width = 0.2,
    position = position_dodge(width = 0.5)
  ) +
  geom_hline(yintercept = 1, color = "dodgerblue", linetype = "dashed", linewidth = 1) +
  facet_wrap(~NA_L1NAME) +
  theme_bw() +
  ylab("Value") +
  xlab("Diversity Index") +
  scale_x_discrete(labels = c(
    "beta_S" = "βS",
    "beta_S_PIE" = "βSPIE",
    "beta_S_C" = "βC"
  )) +
  scale_color_brewer(palette = "Set1") +
  theme(legend.position = "right")








######### MODEL ##########
## the model

library(lme4)
library(lmerTest)
library(emmeans)
library(ggplot2)
library(dplyr)
library(glmmTMB)
library(broom.mixed)   # tidy for mixed models
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(flextable)
library(officer)
library(knitr)
library(kableExtra)

for (idx in unique(beta_results$index)) {
  cat("\n===== Index:", idx, "=====\n")
  
  df <- beta_results %>% filter(index == idx)
  
  # Mixed model (YEAR as random intercept)
  glmer_model <- glmmTMB(value ~ NA_L1NAME + SEASON + (1 | YEAR), 
                         family = Gamma(link = "log"), 
                         data = df)
  print(check_model(glmer_model))
  print(summary(glmer_model))
  
  # Tukey post-hoc for NA_L1NAME
  emm <- emmeans(glmer_model, ~ NA_L1NAME)
  
  tukey_pairs <- as.data.frame(pairs(emm, adjust = "tukey")) %>%
    select(contrast, estimate, p.value)
  
  tukey_ci <- as.data.frame(confint(pairs(emm, adjust = "tukey"))) %>%
    rename(
      lower.CL = asymp.LCL,
      upper.CL = asymp.UCL
    ) %>%
    select(contrast, lower.CL, upper.CL)
  
  
  
  # Join p-values and CI into one clean table
  tukey_df <- left_join(tukey_pairs, tukey_ci, by = "contrast") %>%
    rename(
      Comparison = contrast,
      diff = estimate,
      lwr = lower.CL,
      upr = upper.CL,
      p_value = p.value
    )
  
  cat("\n--- Tukey post-hoc for NA_L1NAME (with CI and p-values) ---\n")
  print(tukey_df)
  
  # Tukey plot
  p <- ggplot(tukey_df, aes(x = reorder(Comparison, diff), y = diff)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    coord_flip() +
    labs(
      title = paste("Tukey HSD –", idx),
      x = "Comparison",
      y = "Difference in Marginal Means"
    ) +
    theme_minimal(base_size = 14)
  
  print(p)
  
  #Tidy fixed effects from your lmer model
  tbl <- broom.mixed::tidy(glmer_model, effects = "fixed") %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))
  
  # 2. Convert to table and export to Word
  ft <- flextable(tbl) 
  doc <- read_docx() %>% body_add_flextable(ft)
  print(doc, target = paste0("model_fixed_effects",idx,".docx"))
  
}



### Supplementary analysis with classic beta metrics
### for just the summer period

calculate_beta_ecoregion_sf <- function(
    wet_long,
    env_eco,
    n_reps = 999,
    C_target_gamma = 0.75,
    max_anchor_attempts = 1000,
    include_mob = TRUE,
    include_classic = TRUE
) {
  
  if (include_classic && !requireNamespace("vegan", quietly = TRUE)) {
    stop("Package 'vegan' is required when include_classic = TRUE.")
  }
  
  calc_classic_beta <- function(comm_mat, eco, anchor_year, anchor_id, rep_id) {
    if (nrow(comm_mat) < 2) {
      return(NULL)
    }
    
    bray_vals <- as.numeric(vegan::vegdist(comm_mat, method = "bray"))
    jaccard_vals <- as.numeric(vegan::vegdist(comm_mat, method = "jaccard", binary = TRUE))
    sorensen_vals <- as.numeric(vegan::vegdist(comm_mat, method = "bray", binary = TRUE))
    
    tibble(
      replicate = rep_id,
      NA_L1NAME = eco,
      YEAR = anchor_year,
      ANCHOR_ID = anchor_id,
      scale = "beta_classic",
      index = c("Bray-Curtis", "Jaccard", "Sorensen"),
      sample_size = NA_real_,
      effort = 5,
      gamma_coverage = NA_real_,
      value = c(
        mean(bray_vals, na.rm = TRUE),
        mean(jaccard_vals, na.rm = TRUE),
        mean(sorensen_vals, na.rm = TRUE)
      )
    )
  }
  
  results_df <- tibble()
  
  env_eco <- env_eco %>%
    filter(LOCALITY_ID %in% wet_long$LOCALITY_ID) %>%
    filter(!is.na(LATITUDE), !is.na(LONGITUDE)) %>%
    filter(
      LATITUDE >= -90, LATITUDE <= 90,
      LONGITUDE >= -180, LONGITUDE <= 180
    )
  
  ecoregions <- unique(env_eco$NA_L1NAME)
  
  for (eco in ecoregions) {
    message("Processing ecoregion: ", eco)
    
    env_sub <- env_eco %>% filter(NA_L1NAME == eco)
    if (nrow(env_sub) < 5) next
    
    env_sub_sf <- st_as_sf(env_sub, coords = c("LONGITUDE", "LATITUDE"), crs = 4326) %>%
      st_transform(crs = 5070)
    
    dist_matrix <- st_distance(env_sub_sf)
    mean_dist_num <- mean(as.numeric(dist_matrix[upper.tri(dist_matrix)]), na.rm = TRUE)
    
    replicate_counter <- 0
    
    while (replicate_counter < n_reps) {
      anchor_attempts <- 0
      valid <- FALSE
      
      while (!valid && anchor_attempts < max_anchor_attempts) {
        anchor_idx <- sample(seq_len(nrow(env_sub_sf)), 1)
        anchor_id <- env_sub$LOCALITY_ID[anchor_idx]
        
        valid_years <- wet_long %>%
          filter(LOCALITY_ID == anchor_id) %>%
          pull(YEAR) %>%
          unique()
        
        if (length(valid_years) == 0) {
          anchor_attempts <- anchor_attempts + 1
          next
        }
        
        anchor_year <- sample(valid_years, 1)
        
        dists <- st_distance(env_sub_sf[anchor_idx, ], env_sub_sf)
        dists_num <- as.numeric(dists)
        candidates <- which(dists_num <= mean_dist_num & dists_num > 0)
        candidates <- unlist(candidates, use.names = FALSE)
        
        candidate_idxs_filtered <- candidates[
          vapply(candidates, function(i) {
            lid <- env_sub$LOCALITY_ID[i]
            n_events <- wet_long %>%
              filter(LOCALITY_ID == lid, YEAR == anchor_year) %>%
              distinct(SAMPLING_EVENT_IDENTIFIER) %>%
              nrow()
            n_events >= 5
          }, logical(1))
        ]
        
        if (length(candidate_idxs_filtered) < 4) {
          anchor_attempts <- anchor_attempts + 1
          next
        }
        
        other_idxs <- sample(candidate_idxs_filtered, 4)
        selected_idxs <- c(anchor_idx, other_idxs)
        selected_ids <- env_sub$LOCALITY_ID[selected_idxs]
        
        events_by_loc <- purrr::map(selected_ids, function(lid) {
          x <- wet_long %>%
            filter(LOCALITY_ID == lid, YEAR == anchor_year) %>%
            distinct(SAMPLING_EVENT_IDENTIFIER) %>%
            pull(SAMPLING_EVENT_IDENTIFIER)
          
          if (length(x) < 5) return(NULL)
          sample(x, size = 5)
        })
        
        if (any(sapply(events_by_loc, is.null))) {
          anchor_attempts <- anchor_attempts + 1
          next
        }
        
        events_map <- tibble(
          LOCALITY_ID = rep(selected_ids, each = 5),
          SAMPLING_EVENT_IDENTIFIER = unlist(events_by_loc)
        )
        
        comm_matrix_df <- wet_long %>%
          filter(
            YEAR == anchor_year,
            SAMPLING_EVENT_IDENTIFIER %in% events_map$SAMPLING_EVENT_IDENTIFIER,
            LOCALITY_ID %in% selected_ids
          ) %>%
          group_by(LOCALITY_ID, COMMON_NAME) %>%
          summarise(abundance = sum(OBSERVATION_COUNT, na.rm = TRUE), .groups = "drop") %>%
          pivot_wider(
            names_from = COMMON_NAME,
            values_from = abundance,
            values_fill = list(abundance = 0)
          )
        
        if (!"LOCALITY_ID" %in% names(comm_matrix_df)) {
          anchor_attempts <- anchor_attempts + 1
          next
        }
        
        comm_mat <- comm_matrix_df %>%
          column_to_rownames("LOCALITY_ID") %>%
          as.data.frame(stringsAsFactors = FALSE)
        
        if (ncol(comm_mat) > 0) {
          comm_mat <- comm_mat[, colSums(comm_mat) > 0, drop = FALSE]
        }
        
        if (
          nrow(comm_mat) == 5 &&
          ncol(comm_mat) > 0 &&
          all(rowSums(comm_mat) > 0)
        ) {
          valid <- TRUE
        } else {
          anchor_attempts <- anchor_attempts + 1
        }
      }
      
      if (!valid) break
      
      out <- list()
      
      if (include_mob) {
        comm_div <- tryCatch({
          calc_comm_div(
            comm_mat,
            index = c("S", "S_PIE", "S_C"),
            extrapolate = TRUE,
            scales = c("beta"),
            C_target_gamma = C_target_gamma
          )
        }, error = function(e) NULL)
        
        if (!is.null(comm_div) && nrow(comm_div) > 0) {
          out[["mob"]] <- as_tibble(comm_div) %>%
            mutate(
              replicate = replicate_counter + 1,
              NA_L1NAME = eco,
              YEAR = anchor_year,
              ANCHOR_ID = anchor_id
            )
        }
      }
      
      if (include_classic) {
        out[["classic"]] <- calc_classic_beta(
          comm_mat = comm_mat,
          eco = eco,
          anchor_year = anchor_year,
          anchor_id = anchor_id,
          rep_id = replicate_counter + 1
        )
      }
      
      out_df <- purrr::list_rbind(out)
      
      if (is.null(out_df) || nrow(out_df) == 0) {
        next
      }
      
      results_df <- bind_rows(results_df, out_df)
      replicate_counter <- replicate_counter + 1
    }
    
    if (replicate_counter < n_reps) {
      warning(sprintf("Only %d valid replicates for ecoregion %s", replicate_counter, eco))
    }
  }
  
  results_df
}

beta_results_summer <- calculate_beta_ecoregion_sf(
  wet_eco_summer,
  env_eco,
  n_reps = 99,
  max_anchor_attempts = 99
)

index_levels <- c("beta_S", "beta_S_PIE", "beta_S_C", "Bray-Curtis", "Jaccard", "Sorensen")


beta_results_summer2 <- beta_results_summer %>%
  mutate(
    index = factor(index, levels = index_levels),
    scale = factor(scale, levels = c("beta", "beta_classic")),
    eco_scale = paste(NA_L1NAME, ifelse(scale == "beta", "MoB", "Classic"), sep = " | ")
  )

wet_div_error <- beta_results_summer2 %>%
  group_by(NA_L1NAME, scale, eco_scale, index) %>%
  summarize(
    mean = mean(value, na.rm = TRUE),
    lower = quantile(value, 0.025, na.rm = TRUE),
    upper = quantile(value, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(beta_results_summer2, aes(x = index, y = value, color = NA_L1NAME)) +
  geom_jitter(alpha = 0.25, width = 0.15) +
  geom_point(
    data = wet_div_error,
    aes(x = index, y = mean),
    inherit.aes = FALSE,
    size = 2.5,
    shape = 21,
    fill = "darkred"
  ) +
  geom_errorbar(
    data = wet_div_error,
    aes(x = index, ymin = lower, ymax = upper),
    inherit.aes = FALSE,
    color = "darkred",
    linewidth = 1,
    width = 0.2
  ) +
  geom_hline(
    data = data.frame(scale = "beta"),
    aes(yintercept = 1),
    inherit.aes = FALSE,
    color = "dodgerblue",
    linetype = "dashed",
    linewidth = 0.9
  ) +
  facet_wrap(
    NA_L1NAME ~ scale,
    scales = "free",
    labeller = labeller(
      NA_L1NAME = label_wrap_gen(width = 18),
      scale = as_labeller(c(
        beta = "MoB β-Diversity",
        beta_classic = "Classic β-Diversity"
      ))
    )
  )+
  xlab("Beta Diversity Index") +
  ylab("Value") +
  scale_x_discrete(labels = c(
    "beta_S" = "βS",
    "beta_S_PIE" = "βSPIE",
    "beta_S_C" = "βC",
    "Bray-Curtis" = "Bray-Curtis",
    "Jaccard" = "Jaccard's Index",
    "Sorensen" = "Sorensen's Index"
  )) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 6),
    panel.spacing = unit(0.8, "lines"),
    legend.position = "none"
  )




