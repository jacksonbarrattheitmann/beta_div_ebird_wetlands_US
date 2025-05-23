# Generating single site-species matrix for all 207 wetlands

library(dplyr)
library(tidyverse)
library(vegan)
library(mobr)

# First let's read in all the GIW data as a single dataframe

wet_dat <- readRDS("Data/eBird_local_alpha_level/ebird_alpha_wetlands_raw.RDS")

# The ENV dataset

env <- readRDS("Data/earth_engine_env_data/env_matrix.RDS")

wet_coords <- wet_dat %>%
  select(LOCALITY_ID, LONGITUDE, LATITUDE) %>%
  distinct()

# let's subset to make this managebale

env_filt <- env %>%
  select(LOCALITY_ID, NA_L1NAME)

# Combine wet_coords and env_filt
env_filt <- env_filt %>%
  inner_join(wet_coords, by = "LOCALITY_ID")


# let's create a site-species matrix for all 207 GIWs from the minimum number of 
# checklists that we have from each GIW

# First, let's check and see how many checklists we have at each GIW

sum_check <- wet_dat %>%
  group_by(LOCALITY_ID) %>%
  summarize(num_check = length(unique(SAMPLING_EVENT_IDENTIFIER)))

# Max = 3550
# Min = 73



# We have a few options
# 1) grab only 73 checklists from all the GIWs, and collapse them all, calculate beta 1 time
# 2) Bootstrapping approach, randomly sub-sample 25-50 and calculate beta each time

# Let's just try #1 for now, and use mobr to analyze

# First thing is we need to create a sample of SAMPLING_EVENT_IDENTIFIERS
# with length 73, for each LOCALITY_ID (eBird hotspot)

check_samples <- wet_dat %>%
  select(LOCALITY_ID, SAMPLING_EVENT_IDENTIFIER) %>%
  distinct(LOCALITY_ID, SAMPLING_EVENT_IDENTIFIER) %>%
  group_by(LOCALITY_ID) %>%
  slice_sample(n = 73)

# check to make sure this worked
# if this is an empty data frame we are good to go
test <- check_samples %>%
  group_by(LOCALITY_ID) %>%
  summarize(num = length(unique(SAMPLING_EVENT_IDENTIFIER))) %>%
  filter(num != 73)

# Now we can filter our wet_dat to only include SAMPLING_EVENT_IDENTIFIERS in
# check_smaples

wet_dat_73 <- wet_dat %>%
  filter(SAMPLING_EVENT_IDENTIFIER %in% check_samples$SAMPLING_EVENT_IDENTIFIER) %>%
  group_by(LOCALITY_ID, COMMON_NAME) %>%
  summarize(
    count = sum(OBSERVATION_COUNT)
  ) %>%
 pivot_wider(names_from = COMMON_NAME, values_from = count, values_fill = 0)

# getting rid of LOCALITY_ID as a column, just species
wet_comm <- wet_dat_73[ , 2:514]

# Have to make sure the row.names are equivlaent
row.names(wet_comm) <- wet_dat_73$LOCALITY_ID


# WE have 1 ecoregion with only 1 site, likley need to filter it out for mobr
# L879018 is the LOCALITY_ID

wet_comm_filt <- wet_dat_73 %>%
  filter(LOCALITY_ID != "L879018")

wet_comm_filt <- wet_comm_filt %>%
  column_to_rownames(var = "LOCALITY_ID")

# Make the row names the same for env_filt and wet_coords
env_filt <- env_filt %>%
  filter(LOCALITY_ID != "L879018")

env_filt <- env_filt %>%
  column_to_rownames(var = "LOCALITY_ID")

## Try just making the calc_biodiv object
## for exploration of the data

biodiv_comm_alpha <- calc_biodiv(wet_comm, groups = env$LOCALITY_ID, 
                                 index = c('N', 'S', 'S_n', 'S_PIE'), 
                                 effort = 50,
                                 extrapolate = FALSE)

# Need to join this object with the env data, to plot the S, S_n, and S_PIE
# within each ecoregion

div_env <- biodiv_comm_alpha %>%
  rename(LOCALITY_ID = group) %>%
  inner_join(env, by = "LOCALITY_ID")

# Calculating the means for plotting
mean_data <- div_env %>%
  group_by(index, NA_L1NAME) %>%
  summarise(mean_value = mean(value), .groups = "drop")

# Now we can plot each of the metrics within the ecoregion, along with the mean
# so that we can try and interpret the betas below

ggplot() + geom_point(data = div_env, aes(x = NA_L1NAME, y = value, color = NA_L1NAME)) +
  geom_point(data = mean_data, aes(x = NA_L1NAME, y = mean_value), color = "red", size = 3) +
  facet_wrap(~index, scales = "free_y") +
  theme(axis.text = element_text(angle = 90))


## Trying to use Vegan to calculate beta diversity patterns among all 
## GIWs then maybe we can plot these against each other by ecoregion

whit_beta <- betadiver(wet_comm_filt, "w")

beta_NA <- betadisper(whit_beta, env_filt$NA_L1NAME, type = "median")

mod <- anova(beta_NA)

permutest(beta_NA, pairwise = TRUE, permutations = 99)

mod.HSD <- TukeyHSD(beta_NA)

plot(mod.HSD)

plot(beta_NA)

plot(beta_NA, ellipse = TRUE, hull = FALSE, conf = 0.90) # 90% data ellipse

boxplot(beta_NA)


#### THIS MAY BE UNNECASSARY and NOT USEFUL for BETA Diversity comaprions among GIWs
# Now create the mob object
wet_mob <- make_mob_in(wet_comm_filt, plot_attr = env_filt, coord_names = c('LONGITUDE', 'LATITUDE'))

# Take a long time to run, so I have this saved as an RDS below for loading into R

# beta_results <- get_delta_stats(wet_mob, env_var = "NA_L1NAME", stats = "betas", type = 'discrete',
#                                 log_scale = TRUE)

# saveRDS(beta_results, "Data/beta_results.RDS")

beta_results <- readRDS("Data/beta_results.RDS")

plot(beta_results, stat = 'b1')

betas <- beta_results$S_df



