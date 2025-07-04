# analysis script

library(dplyr)
library(tidyverse)
library(vegan)
library(mobr)
library(ggplot2)
library(betapart)
library(lme4)
library(performance)

######### OBJETICVE 1 - CONTINENTAL SCALE ########
######### ALL GIWs beta ##############

wet_all <- readRDS("Intermediate_data/wet_comm_all_summarized.RDS") %>%
  column_to_rownames(var = "LOCALITY_ID")

betas_ALL <- calc_comm_div(wet_all, index = c('S','S_n','S_PIE', 'S_C'),
                             extrapolate = TRUE, effort = 25, scales = c("beta"), C_target_gamma = 0.75)

# Classic Whittaker's beta diversity (global)
sitespp_pa <- (wet_all > 0) * 1
alpha <- rowSums(sitespp_pa)
mean_alpha <- mean(alpha)
gamma <- sum(colSums(sitespp_pa) > 0)

beta_whittaker <- as.data.frame(gamma / mean_alpha - 1)

# Jaccard's
beta_dist <- vegdist(sitespp_pa, method = "jaccard", binary = TRUE)
mean(beta_dist)

beta_jac <- as.data.frame(mean(beta_dist))

# Whit and Jac toghther in a single df
beta_traditional <- cbind(beta_jac, beta_whittaker) %>%
  rename("Whitaker_beta" = "gamma/mean_alpha - 1",
         "Jaccard's Index" = "mean(beta_dist)")

#combining them all toghther
betas_fin <- cbind(betas, beta_traditional) %>%
  pivot_longer(cols = c("beta_S", "beta_S_n", "beta_S_PIE", "beta_S_PIE", "beta_S_C", "Whitaker_beta", "Jaccard's Index"), 
               names_to = "index", values_to = "value")
########## PLOTS for OBJECTIVE 1 ##################

ggplot(data = betas_ALL) +
  geom_col(aes(x = index, y = value), fill = "lightblue") + theme_bw() +
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


######### OBJECTIVE 3 - GIW scale ############
######### Within GIW beta diversity analysis ##########

#read in env and wet_comm data
env <- readRDS("Intermediate_data/env_all_wet_check.RDS")

wet_dat <- readRDS("Intermediate_data/all_wet_check_summarized.RDS")


## make SAMPLING_EVENT_IDENTIFIER into the row names

env <- column_to_rownames(env, var = "SAMPLING_EVENT_IDENTIFIER")

# do the same for wet_dat
wet_dat <- column_to_rownames(wet_dat, var = "SAMPLING_EVENT_IDENTIFIER")

# need to remove LOCALITY_ID column
wet_dat <- wet_dat[, 2:529]

######### now time to run mobr ##################
calc_C_target(wet_dat)

wet_div <- tibble(wet_dat) %>%
  group_by(group = env$LOCALITY_ID) %>%
  group_modify(~ calc_comm_div(.x,
                               index = c('S', 'S_n','S_PIE', 'S_C'),
                               extrapolate = TRUE,
                               effort = 25,
                               scales = "beta",
                               C_target_gamma = 0.75))


wet_div_error <- wet_div %>%
  group_by(index) %>%
  summarize(
    mean = mean(value, na.rm = TRUE),
    lower = quantile(value, 0.025, na.rm = TRUE),
    upper = quantile(value, 0.975, na.rm = TRUE)
  )


# Plotting the results
ggplot() +
  geom_jitter(data = wet_div, aes(x = index, y = value), alpha = 0.25, width = 0.2) +
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


# Plotting as geom_col
ggplot(wet_div_error, aes(x = index, y = mean)) +
  geom_col(fill = "skyblue") +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3, color = "red") +
  geom_hline(yintercept = 1, color = "darkred", linetype = "dashed", linewidth = 1) +
  theme_bw() +
  ylab("Value") +
  xlab("Diversity Index") +
  scale_x_discrete(labels = c(
    "beta_S" = "βS",
    "beta_S_n" = "βSn",
    "beta_S_PIE" = "βSPIE", 
    "beta_S_C" = "βC"))


# need to add the GIW area per state as an explanatory variable

###### building a model to explain variation in aggregation ############
wet_div_betas_filt <- wet_div_betas %>%
  rename(LOCALITY_ID = group)

# first I need to append the beta stats into a table with the env data
data_mod <- env %>%
  distinct(LOCALITY_ID, .keep_all = TRUE) %>%
  full_join(wet_div_betas_filt, by = "LOCALITY_ID")


mod <- lm(value ~ built_wet + flooded_vegetation_wet + area_sqkm + shan_wet + evi_mean + water_25km, data = data_mod)
summary(mod)
check_model(mod)

# wetland area
ggplot(data = data_mod) +
  geom_point(aes(x = log(area_sqkm*1000), y = value)) +
  geom_smooth(aes(x = log(area_sqkm*1000), y = value), method = "lm") +
  theme_bw()

# EVI
ggplot(data = data_mod) +
  geom_point(aes(x = evi_mean, y = value)) +
  geom_smooth(aes(x = evi_mean, y = value), method = "lm") +
  theme_bw()

# habitat heterogeneity
ggplot(data = data_mod) +
  geom_point(aes(x = shan_wet, y = value)) +
  geom_smooth(aes(x = shan_wet, y = value), method = "lm") +
  theme_bw()



