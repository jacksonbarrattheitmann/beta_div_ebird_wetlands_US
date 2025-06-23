# Generating single site-species matrix for all 207 wetlands

library(dplyr)
library(tidyverse)
library(vegan)
library(mobr)
library(ggplot2)
library(betapart)

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

# Have to make sure the row.names are equivalent
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

## Let's use the mobr framework as Dan suggested
## 

wet_mob <- make_mob_in(wet_comm_filt, env_filt, coord_names = c("LONGITUDE", "LATITUDE"))
calc_beta_div(wet_comm_filt, 'S_C')  # this is beta_C for the entire matrix
calc_beta_div(wet_comm_filt, 'S_C', C_target_gamma = 0.5)  # at 50% coverage 

# need to use calc_comm_div


wet_div <- tibble(wet_comm_filt) %>% 
  group_by(group = env_filt$NA_L1NAME) %>% 
  group_modify(~ calc_comm_div(.x, index = c('N','S','S_n', 'S_PIE', 'S_C'), effort = 25,
                               extrapolate = TRUE), scales = c("alpha", "gamma", "beta"))

# then plot using plot_comm_div
plot_comm_div(wet_div, multi_panel = FALSE)


## Trying to make the plot in ggplot to customize visuals

wet_div %>%
  filter(index == "S_C" | index == "beta_S_C") %>%
ggplot() +
  geom_boxplot(aes(x = group, y = value, color = group)) +
  facet_wrap(~scale, scales = 'free_y') +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        legend.position = "none")


## printing the beta_S_C data

wet_div %>%
  filter(index == "beta_S_C")



# Eastern Temp Forests first

env_filt_ETF <- env %>%
  group_by(NA_L1NAME) %>%
  filter(NA_L1NAME == "EASTERN TEMPERATE FORESTS")

wet_comm_ETF <- wet_dat_73 %>%
  filter(LOCALITY_ID %in% env_filt_ETF$LOCALITY_ID)

env_filt_ETF <- env_filt_ETF %>%
  column_to_rownames(var = "LOCALITY_ID")

wet_comm_ETF <- wet_comm_ETF %>%
  column_to_rownames(var = "LOCALITY_ID")

# Compute gamma (522)
gamma_ETF <- as.numeric(length(wet_comm_ETF))

# Compute alpha bar (110)
alpha_bar <- wet_comm_ETF %>%
  mutate(richness = rowSums(. > 0)) %>% 
  summarise(mean_richness = mean(richness)) %>%
  pull(mean_richness) %>%
  as.numeric()
  
# Now we can try to compute Whitaker's Beta
# for just EASTERN TEMPERATE FORESTS as a test

beta_ETF <- betadiver(wet_comm_ETF, "w")
# mean beta between all sites
# 0.37 = low turnover
median(beta_ETF)

# Try Sorenson's
beta_ETF_S <- betadiver(wet_comm_ETF, method = 11)

mean(beta_ETF_S)

# breaking up nestedness and turnover
bas_S <- nestedbetajac(wet_comm_ETF)
bas_S

nestedchecker(wet_comm_ETF)

# Here we can measure what the mean distance is (Bray-Curtis)
dist_beta <- as.data.frame(mean(vegdist(wet_comm_ETF, "bray")))

# Then compute a null model to compare
meandist <- function(x) mean(vegdist(x, "bray"))
mbc1 <- oecosimu(wet_comm_ETF, meandist, "r2dtable")

# PLotting
sim_vals <- as.numeric(mbc1$oecosimu$simulated)
obs_val <- mbc1$statistic
null_mean <- mean(sim_vals)

# Create histogram with breaks that include observed value
breaks_seq <- seq(min(sim_vals), max(c(sim_vals, obs_val)) + 0.05, length.out = 20)

hist(sim_vals,
     main = "Null distribution of mean Bray-Curtis dissimilarity",
     xlab = "Simulated mean dissimilarity",
     col = "lightgray",
     border = "white",
     breaks = breaks_seq,
     xlim = c(min(breaks_seq), max(breaks_seq)))

# Add vertical lines
abline(v = obs_val, col = "red", lwd = 2)
abline(v = null_mean, col = "blue", lwd = 2, lty = 2)

# Add legend with both values
legend("top", 
       legend = c(paste("Observed =", round(obs_val, 4)),
                  paste("Null mean =", round(null_mean, 4))),
       col = c("red", "blue"), 
       lwd = 2, 
       lty = c(1, 2),
       bty = "n")


# now we can try betadispersion in vegan

beta_ETF_dis <- betadisper(beta_ETF, env_filt_ETF$NA_L2NAME, "median")

plot(beta_ETF_dis)

anova(beta_ETF_dis)
TukeyHSD((beta_ETF_dis))












## Not sure how useful all of this is
beta_NA <- betadisper(whit_beta, env_filt$NA_L1NAME, type = "median")

mod <- anova(beta_NA)

permutest(beta_NA, pairwise = TRUE, permutations = 99)

mod.HSD <- TukeyHSD(beta_NA)

plot(mod.HSD)

plot(beta_NA)

plot(beta_NA, ellipse = TRUE, hull = FALSE, conf = 0.90) # 90% data ellipse

boxplot(beta_NA)





