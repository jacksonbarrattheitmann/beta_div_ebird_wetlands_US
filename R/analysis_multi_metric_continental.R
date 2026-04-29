# analysis script


######### OBJECTIVE 1 - CONTINENTAL SCALE ########
######### ALL GIWs beta per season ##############

wet_all_summer <- readRDS("Intermediate_data/all_wet_check_for_analysis_SUMMER.RDS")

wet_all_fall <- readRDS("Intermediate_data/all_wet_check_for_analysis_FALL.RDS")

wet_all_spring <- readRDS("Intermediate_data/all_wet_check_for_analysis_SPRING.RDS")

wet_all_winter <- readRDS("Intermediate_data/all_wet_check_for_analysis_WINTER.RDS")

env_summer <- readRDS("Data/old_data/earth_engine_env_data/env_matrix.RDS")

# need spatial coords for the LOCALITY_IDs
wet_coords <- readRDS("Intermediate_data/locality_ids_long_lat.RDS")

env <- readRDS("Data/old_data/earth_engine_env_data/env_matrix.RDS")

env <- env %>%
  left_join(wet_coords, by = "LOCALITY_ID")

## Creating a dataframe with just the LOCALITY_ID and ECOREGION
## to append to the wet_all df

ecoregion_loc <- env %>%
  select(LOCALITY_ID, NA_L1NAME)

wet_all_summer <- wet_all_summer %>%
  inner_join(ecoregion_loc, by = "LOCALITY_ID")

wet_all_fall <- wet_all_fall %>%
  inner_join(ecoregion_loc, by = "LOCALITY_ID")

wet_all_spring <- wet_all_spring %>%
  inner_join(ecoregion_loc, by = "LOCALITY_ID")

wet_all_winter <- wet_all_winter %>%
  inner_join(ecoregion_loc, by = "LOCALITY_ID") 
 

##### FUNCTION for calculating beta across ALL SITES (ONLY FOR SUPPLEMENT)

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
        index = c("S", "S_PIE", "S_C"),
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


betas_ALL_summer <- calc_beta_by_year(wet_all_summer, effort = 5, seed = 5)

betas_ALL_summer <-  betas_ALL_summer %>%
  mutate(SEASON = "Summer")

betas_ALL_winter <- calc_beta_by_year(wet_all_winter, effort = 5, seed = 5)

betas_ALL_winter <-  betas_ALL_winter %>%
  mutate(SEASON = "Winter")

betas_ALL_spring <- calc_beta_by_year(wet_all_spring, effort = 5, seed = 5)

betas_ALL_spring <- betas_ALL_spring %>%
  mutate(SEASON = "Spring")

betas_ALL_fall <- calc_beta_by_year(wet_all_fall, effort = 5, seed = 5)

betas_ALL_fall <- betas_ALL_fall %>%
  mutate(SEASON = "Fall")

betas_ALL <- rbind(betas_ALL_winter, betas_ALL_spring, betas_ALL_summer, betas_ALL_fall)

# SUpplement plot

betas_boot %>%
  filter(SEASON == "Summer") %>%
ggplot() +
  geom_point(aes(x = factor(index, levels = c("beta_S", "beta_S_PIE", "beta_S_C")), y = value, color = SEASON)) + theme_bw() +
  geom_hline(yintercept = 1, color = "dodgerblue", linetype = "dashed", linewidth = 1) +
  xlab("Beta Diversity Index") +
  ylab("Value") +
  facet_wrap(~SEASON) +
  scale_x_discrete(  labels = c(
    "beta_S" = "βS",
    "beta_S_PIE" = "βSPIE", 
    "beta_S_C" = "βC"
  )) +
  theme(legend.position = "right",
        axis.title.x = element_text(
          margin = margin(t = 15)),
        axis.title.y = element_text(
          margin = margin(r = 15))) +
  scale_color_brewer(palette = "Set1")


ggsave("Fig1_betas_CONTINENTAL_SCALE.png", width = 6, height = 4,
       bg = "transparent")


## OBJECTIVE CONTINENTAL SCALE - with SENSITIVITY FOR SITES per ECOREGION #############
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
          index = c("S", "S_PIE", "S_C"),
          extrapolate = TRUE,
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

betas_boot_summer <- calc_beta_by_year_boot(wet_all_summer, effort = 5, n_boot = 99, seed = NULL)

betas_boot_summer <- betas_boot_summer %>%
  mutate(SEASON = "Summer")

betas_boot_fall <- calc_beta_by_year_boot(wet_all_fall, effort = 5, n_boot = 99, seed = NULL)

betas_boot_fall <- betas_boot_fall %>%
  mutate(SEASON = "Fall")

betas_boot_spring <- calc_beta_by_year_boot(wet_all_spring, effort = 5, n_boot = 99, seed = NULL)

betas_boot_spring <- betas_boot_spring %>%
  mutate(SEASON = "Spring")

betas_boot_winter <- calc_beta_by_year_boot(wet_all_winter, effort = 5, n_boot = 99, seed = NULL)

betas_boot_winter <- betas_boot_winter %>%
  mutate(SEASON = "Winter")

## Let's rbind all of these dfs toghther for plotting

betas_boot <- rbind(betas_boot_winter, betas_boot_spring, betas_boot_summer, betas_boot_fall)

#saveRDS(betas_boot, "Intermediate_data/beta_results_continental.RDS")

betas_boot <- readRDS("Intermediate_data/beta_results_continental.RDS")

betas_sum <- betas_boot %>%
  filter(SEASON == "Summer")

wet_div_error <- betas_sum %>%
  group_by(index) %>%
  summarize(
    mean = mean(value, na.rm = TRUE),
    lower = quantile(value, 0.025, na.rm = TRUE),
    upper = quantile(value, 0.975, na.rm = TRUE)
  )

# Plot for the main text
ggplot(data = betas_sum) +
  geom_jitter(aes(x = factor(index, levels = c("beta_S", "beta_S_PIE", "beta_S_C")), y = value), alpha = 0.25) +
  geom_point(data = wet_div_error,aes(x = factor(index, levels = c("beta_S", "beta_S_PIE", "beta_S_C")),
      y = mean
    ),
    size = 3,
    shape = 21,
    fill = "darkred"
  ) +
  geom_errorbar (data = wet_div_error,aes(x = factor(index, levels = c("beta_S", "beta_S_PIE", "beta_S_C")),
      ymin = lower,
      ymax = upper
    ),
    color = "darkred",
    linewidth = 1.2,
    width = 0.2
  ) +
  geom_hline(yintercept = 1, color = "dodgerblue", linetype = "dashed", linewidth = 1) +
  xlab("Beta Diversity Index") +
  ylab("Value") +
  scale_x_discrete(  labels = c(
    "beta_S" = "βS",
    "beta_S_PIE" = "βSPIE", 
    "beta_S_C" = "βC"
  )) +
  theme(legend.position = "right",
        axis.title.x = element_text(
          margin = margin(t = 15)),
        axis.title.y = element_text(
          margin = margin(r = 15))) +
    theme_bw()

## Model results for SEASON

for (idx in unique(betas_boot$index)) {
  cat("\n===== Index:", idx, "=====\n")
  
  df <- betas_boot %>% filter(index == idx)
  
  mod <- glmmTMB(
    value ~ SEASON +
      (1 | YEAR), # random intercept for YEAR#     
    family = Gamma(link = "log"),
    data = df
  )
  print(summary(mod))
  print(check_model(mod))
  
  #Tidy fixed effects from your lmer model
  tbl <- broom.mixed::tidy(mod, effects = "fixed") %>%
    mutate(across(where(is.numeric), ~ round(.x, 3)))
  
  # 2. Convert to table and export to Word
  ft <- flextable(tbl) 
  doc <- read_docx() %>% body_add_flextable(ft)
  print(doc, target = paste0("model_fixed_effects",idx,".docx"))
  
}







# supplement figure
betas_boot %>%
  distinct(SEASON, YEAR, .keep_all = TRUE) %>%
ggplot() +
  geom_point(aes(x = YEAR, y = n_sites, color = SEASON), position = position_dodge(width = 0.1)) +
  scale_color_brewer(palette = "Set1") +
  ylab("Number of GIWs sampled") +
  theme_bw()

# calculating the means and error bars for plotting
wet_div_error <- betas_boot %>%
  group_by(index, SEASON) %>%
  summarize(
    mean = mean(value, na.rm = TRUE),
    lower = quantile(value, 0.025, na.rm = TRUE),
    upper = quantile(value, 0.975, na.rm = TRUE)
  )

# Plotting the results
ggplot() +
  geom_jitter(
    data = betas_boot,
    aes(
      x = factor(index, levels = c("beta_S", "beta_S_PIE", "beta_S_C")),
      y = value,
      color = SEASON
    ),
    alpha = 0.1,
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.5)
  ) +
  geom_point(
    data = wet_div_error,
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
    data = wet_div_error,
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






######## SUPPLEMENTARY ANALYSIS w/ classic Beta Metrics #######

calc_beta_by_year_boot <- function(
    df,
    effort = 5,
    n_boot = 99,
    seed = NULL,
    include_mob = TRUE,
    include_classic = TRUE
) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  if (include_classic && !requireNamespace("vegan", quietly = TRUE)) {
    stop("Package 'vegan' is required when include_classic = TRUE.")
  }
  
  calc_classic_beta <- function(comm, year_val, boot_id, n_sites) {
    if (nrow(comm) < 2) {
      return(NULL)
    }
    
    bray_vals <- as.numeric(vegan::vegdist(comm, method = "bray"))
    jaccard_vals <- as.numeric(vegan::vegdist(comm, method = "jaccard", binary = TRUE))
    sorensen_vals <- as.numeric(vegan::vegdist(comm, method = "bray", binary = TRUE))
    
    tibble::tibble(
      YEAR = year_val,
      bootstrap = boot_id,
      n_sites = n_sites,
      scale = "beta_classic",
      index = c("Bray-Curtis", "Jaccard", "Sorensen"),
      sample_size = NA_real_,
      effort = effort,
      gamma_coverage = NA_real_,
      value = c(
        mean(bray_vals, na.rm = TRUE),
        mean(jaccard_vals, na.rm = TRUE),
        mean(sorensen_vals, na.rm = TRUE)
      )
    )
  }
  
  dplyr::group_by(df, YEAR) %>%
    dplyr::group_split() %>%
    purrr::map_dfr(function(year_dat) {
      year_val <- year_dat$YEAR[1]
      
      sites_ok <- year_dat %>%
        dplyr::distinct(LOCALITY_ID, SAMPLING_EVENT_IDENTIFIER) %>%
        dplyr::count(LOCALITY_ID, name = "n_checklists") %>%
        dplyr::filter(n_checklists >= effort) %>%
        dplyr::pull(LOCALITY_ID)
      
      year_dat <- year_dat %>%
        dplyr::filter(LOCALITY_ID %in% sites_ok)
      
      if (dplyr::n_distinct(year_dat$NA_L1NAME) < 2) {
        return(NULL)
      }
      
      purrr::map_dfr(seq_len(n_boot), function(b) {
        chosen_sites <- year_dat %>%
          dplyr::distinct(NA_L1NAME, LOCALITY_ID) %>%
          dplyr::group_by(NA_L1NAME) %>%
          dplyr::slice_sample(n = 1) %>%
          dplyr::ungroup() %>%
          dplyr::pull(LOCALITY_ID)
        
        chosen_chk <- year_dat %>%
          dplyr::filter(LOCALITY_ID %in% chosen_sites) %>%
          dplyr::distinct(LOCALITY_ID, SAMPLING_EVENT_IDENTIFIER) %>%
          dplyr::group_by(LOCALITY_ID) %>%
          dplyr::slice_sample(n = effort) %>%
          dplyr::ungroup()
        
        sampled <- year_dat %>%
          dplyr::semi_join(
            chosen_chk,
            by = c("LOCALITY_ID", "SAMPLING_EVENT_IDENTIFIER")
          )
        
        comm <- sampled %>%
          dplyr::group_by(LOCALITY_ID, COMMON_NAME) %>%
          dplyr::summarise(
            abundance = sum(OBSERVATION_COUNT, na.rm = TRUE),
            .groups = "drop"
          ) %>%
          tidyr::pivot_wider(
            names_from = COMMON_NAME,
            values_from = abundance,
            values_fill = 0
          ) %>%
          as.data.frame()
        
        rownames(comm) <- comm$LOCALITY_ID
        comm$LOCALITY_ID <- NULL
        
        out <- list()
        
        if (include_mob) {
          beta_out <- calc_comm_div(
            comm,
            index = c("S", "S_PIE", "S_C"),
            extrapolate = TRUE,
            scales = "beta",
            C_target_gamma = 0.75
          )
          
          out[["mob"]] <- tibble::tibble(
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
        }
        
        if (include_classic) {
          out[["classic"]] <- calc_classic_beta(
            comm = comm,
            year_val = year_val,
            boot_id = b,
            n_sites = length(chosen_sites)
          )
        }
        
        purrr::list_rbind(out)
      })
    })
}

betas_boot_summer <- calc_beta_by_year_boot(wet_all_summer, effort = 5, n_boot = 99, seed = NULL)


index_levels <- c("beta_S", "beta_S_PIE", "beta_S_C", "Bray-Curtis", "Jaccard", "Sorensen")

betas_boot_summer2 <- betas_boot_summer %>%
  mutate(
    index = factor(index, levels = index_levels),
    scale = factor(scale, levels = c("beta", "beta_classic"))
  )

wet_div_error <- betas_boot_summer2 %>%
  group_by(scale, index) %>%
  summarize(
    mean = mean(value, na.rm = TRUE),
    lower = quantile(value, 0.025, na.rm = TRUE),
    upper = quantile(value, 0.975, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(betas_boot_summer2, aes(x = index, y = value)) +
  geom_jitter(alpha = 0.25, width = 0.15) +
  geom_point(
    data = wet_div_error,
    aes(x = index, y = mean),
    inherit.aes = FALSE,
    size = 3,
    shape = 21,
    fill = "darkred"
  ) +
  geom_errorbar(
    data = wet_div_error,
    aes(x = index, ymin = lower, ymax = upper),
    inherit.aes = FALSE,
    color = "darkred",
    linewidth = 1.2,
    width = 0.2
  ) +
  geom_hline(
    data = data.frame(scale = "beta"),
    aes(yintercept = 1),
    inherit.aes = FALSE,
    color = "dodgerblue",
    linetype = "dashed",
    linewidth = 1
  ) +
  facet_wrap(
    ~scale,
    scales = "free",
    labeller = as_labeller(c(
      beta = "MoB β-Diversity",
      beta_classic = "Classic β-Diversity"
    ))
  ) +
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
  theme_bw()
