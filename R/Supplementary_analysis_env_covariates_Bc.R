#### Supplementary analysis for env covariate impacts ####

env_summer_join <- env_summer %>%
  select(LOCALITY_ID, evi_mean, hum_pop_mean, water_25km, flooded_vegetation_25km)

wet_all_summer_new <- wet_all_summer %>%
  left_join(env_summer_join, by = "LOCALITY_ID")
  



calc_beta_by_year_boot <- function(
    df,
    effort = 5,
    n_boot = 99,
    seed = NULL,
    env_vars = c("evi_mean", "hum_pop_mean", "water_25km", "flooded_vegetation_25km"),
    lon_var = "LONGITUDE",
    lat_var = "LATITUDE"
) {
  if (!is.null(seed)) set.seed(seed)
  
  df %>%
    dplyr::group_by(YEAR) %>%
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
      
      if (dplyr::n_distinct(year_dat$NA_L1NAME) < 2) return(NULL)
      
      site_env <- year_dat %>%
        dplyr::distinct(LOCALITY_ID, NA_L1NAME, dplyr::across(dplyr::all_of(c(env_vars, lon_var, lat_var))))
      
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
          dplyr::semi_join(chosen_chk, by = c("LOCALITY_ID", "SAMPLING_EVENT_IDENTIFIER"))
        
        comm <- sampled %>%
          dplyr::group_by(LOCALITY_ID, COMMON_NAME) %>%
          dplyr::summarise(abundance = sum(OBSERVATION_COUNT, na.rm = TRUE), .groups = "drop") %>%
          tidyr::pivot_wider(names_from = COMMON_NAME, values_from = abundance, values_fill = 0) %>%
          as.data.frame()
        
        rownames(comm) <- comm$LOCALITY_ID
        comm$LOCALITY_ID <- NULL
        
        beta_out <- calc_comm_div(
          comm,
          index = c("S", "S_PIE", "S_C"),
          extrapolate = TRUE,
          scales = "beta",
          C_target_gamma = 0.75
        )
        
        env_boot <- site_env %>%
          dplyr::filter(LOCALITY_ID %in% chosen_sites)
        
        env_summary <- purrr::map_dfc(env_vars, function(v) {
          x <- env_boot[[v]]
          tibble::tibble(
            !!paste0("mean_", v) := mean(x, na.rm = TRUE),
            !!paste0("sd_", v) := stats::sd(x, na.rm = TRUE)
          )
        })
        
        geo_summary <- tibble::tibble(mean_geo_dist = NA_real_)
        
        if (nrow(env_boot) > 1) {
          coords <- env_boot %>%
            dplyr::select(dplyr::all_of(c(lon_var, lat_var))) %>%
            as.matrix()
          
          geo_summary$mean_geo_dist <- mean(
            as.numeric(stats::dist(coords)),
            na.rm = TRUE
          )
        }
        
        tibble::tibble(
          YEAR = year_val,
          bootstrap = b,
          n_sites = length(chosen_sites),
          scale = beta_out$scale,
          index = beta_out$index,
          sample_size = beta_out$sample_size,
          effort = beta_out$effort,
          gamma_coverage = beta_out$gamma_coverage,
          value = beta_out$value
        ) %>%
          dplyr::bind_cols(env_summary, geo_summary)
      })
    })
}

beta_c_dat <- calc_beta_by_year_boot(
  wet_all_summer_new,
  effort = 5,
  n_boot = 10,
  seed = 1,
  env_vars = c("evi_mean", "hum_pop_mean", "water_25km", "flooded_vegetation_25km")
) %>%
  dplyr::filter(index == "beta_S_C")

mod_beta_c <- glmmTMB(
  value ~ scale(mean_evi_mean) +
    scale(mean_water_25km) +
    scale(mean_flooded_vegetation_25km) +
    scale(mean_hum_pop_mean)+
    scale(mean_geo_dist) + (1 | YEAR),
  family = Gamma(link = "log"),
  data = beta_c_dat
)

check_model(mod_beta_c)

library(broom.mixed)
library(ggplot2)
library(dplyr)

coef_df <- broom.mixed::tidy(mod_beta_c, effects = "fixed", conf.int = TRUE) %>%
  dplyr::filter(term != "(Intercept)")

ggplot(coef_df, aes(x = estimate, y = reorder(term, estimate))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 3, color = "firebrick3") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2, color = "firebrick3") +
  labs(
    x = "Standardized coefficient",
    y = "Predictor",
    title = expression(beta[C] ~ "environmental drivers")
  ) +
  theme_bw()

