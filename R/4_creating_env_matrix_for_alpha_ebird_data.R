library(dplyr)
library(tidyverse)
library(raster)
library(sf)
library(vegan)

sf::sf_use_s2(FALSE)

# Creating the environmental matrix for eBird data at the alpha level

# First we load in the list of wetlands
wetlands <- read_csv(file = "Data/potential_study_sites_list/potential_wetland_ebird_locations_TRACKING.csv")

# Need to filter out the ones that are not isolated from QA/QC column (we only want Y for Yes)

wetlands <- wetlands %>%
  filter(QAQC == "Y")

wet_type <- wetlands %>%
  dplyr::select(LOCALITY_ID, wet_type)

# Load in all the separate environment csv files, to combine, and clean them of columns we do not need

land_25 <- read_csv(file = "Data/earth_engine_env_data/landcover_est_within_25k_buffer.csv")
land_25 <- dplyr::select(land_25, -"system:index", -".geo", -"label") %>%
  dplyr::rename(LOCALITY_ID = LOCALIT)

land_10 <- read_csv(file = "Data/earth_engine_env_data/landcover_est_within_10k_buffer.csv")
land_10 <- dplyr::select(land_10, -"system:index", -".geo", -"label") %>%
  dplyr::rename(LOCALITY_ID = LOCALIT)

land_5 <- read_csv(file = "Data/earth_engine_env_data/landcover_est_within_5k_buffer.csv")
land_5 <- dplyr::select(land_5, -"system:index", -".geo", -"label") %>%
  dplyr::rename(LOCALITY_ID = LOCALIT)

land_in_wet <- read_csv(file = "Data/earth_engine_env_data/landcover_within_wet.csv")
land_in_wet <- dplyr::select(land_in_wet, -"system:index", -".geo", -"label_wet")

evi_est <- read_csv(file = "Data/earth_engine_env_data/25k_buffer_EVI_estimates.csv")
evi_est <- dplyr::select(evi_est, -"system:index", -".geo")

pop_den <- read_csv(file = "Data/earth_engine_env_data/pop_density_25km_buffer.csv")
pop_den <- dplyr::select(pop_den, -"system:index", -".geo")

# Calculating shannon diversity index of land cover attributes per LOCALITY_ID at the local, and gamma scales
#local scale within wetland site
shan_wet <- data.frame(LOACILITY_ID = land_in_wet[, 1],
                       shan_wet = diversity(land_in_wet[, -1], index = "shannon"))

#gamma scale within 25km of wetland site
shan_gamma_25 <- data.frame(LOACILITY_ID = land_25[, 1],
                         shan_gamma_25 = diversity(land_25[, -1], index = "shannon"))


#gamma scale within 10km of wetland site
shan_gamma_10 <- data.frame(LOACILITY_ID = land_10[, 1],
                         shan_gamma_10 = diversity(land_10[, -1], index = "shannon"))


#gamma scale within 5km of wetland site
shan_gamma_5 <- data.frame(LOACILITY_ID = land_5[, 1],
                         shan_gamma_5 = diversity(land_5[, -1], index = "shannon"))

# Now lets try to calculate the area of each polygon
x <- shapefile('Data/combined_shp_geojson/combined_polygons.shp')
crs(x)
x$area_sqkm <- area(x) / 1000000
area_wet <- as.data.frame(x)
area_wet <- area_wet %>% 
  rename(LOCALITY_ID = LOCALIT)


# Finding the EPA Level II ecoregion for each intersecting polygon
# Loading in the data and making it into an sf object
eco_reg <- read_sf('Data/earth_engine_env_data/epa_level2_eco/NA_CEC_Eco_Level2.shp') %>%
  st_as_sf(crs = st_crs(x))

# Making the projection the same as the wetlands polygon
eco_reg <- st_transform(eco_reg, "+proj=longlat +datum=WGS84 +no_defs")

# making the wetlands polygon, x, into an sf object as well
wet_shp <- x %>%
  st_as_sf()

class(eco_reg)
wet_shp <- st_make_valid(wet_shp)
eco_reg <- st_make_valid(eco_reg)


# bit hacky, but I think it works
pi <- wet_shp %>% 
  st_within(eco_reg) %>%
  as.data.frame() %>%
  left_join(eco_reg %>%
               mutate(col.id=1:nrow(.))) %>%
  dplyr::select(-geometry) %>%
  left_join(wet_shp %>%
              mutate(row.id=1:nrow(.))) %>%
  dplyr::select(-geometry) %>%
  rename(LOCALITY_ID=LOCALIT) %>%
  dplyr::select(LOCALITY_ID, 3:8)


# Now combine them based on unique identifier (LOCALITY_ID)

env_matrix <- evi_est %>%
  full_join(pop_den, by = "LOCALITY_ID") %>%
  left_join(land_in_wet, by = "LOCALITY_ID") %>%
  left_join(land_25, by = "LOCALITY_ID") %>%
  left_join(land_10, by = "LOCALITY_ID") %>%
  left_join(land_5, by = "LOCALITY_ID") %>%
  left_join(area_wet, by = "LOCALITY_ID") %>%
  left_join(pi, by = "LOCALITY_ID") %>%
  left_join(shan_wet, by = "LOCALITY_ID") %>%
  left_join(shan_gamma_25, by = "LOCALITY_ID") %>%
  left_join(shan_gamma_10, by = "LOCALITY_ID") %>%
  left_join(shan_gamma_5, by = "LOCALITY_ID") %>%
  left_join(wet_type, by = "LOCALITY_ID")

# export the env_matrix as a csv

saveRDS(env_matrix, file = "Data/earth_engine_env_data/env_matrix.RDS")




#### Supp table creation


supp <- read.csv("Data/supp_data/Supp_Table.csv")

supp <- supp %>%
  rename("LOCALITY_ID" = "Locality_ID")

supp_2 <- st_as_sf(supp, coords = c("Longitude", "Latitude"), crs = 4326)

#gettign just area from the previosu env matrix
env_areas <- env_matrix[, c(1, 40)]

## adding in US county polygon
us_count <- read_sf('Data/us_counties/cb_2018_us_county_20m.shp') %>%
  st_as_sf(crs = st_crs(x))

us_count <- st_transform(us_count, "+proj=longlat +datum=WGS84 +no_defs")

us_count <- st_make_valid(us_count)

#finding the intersection
pi_2 <- supp_2 %>% 
  st_within(us_count) %>%
  as.data.frame() %>%
  left_join(us_count %>%
              mutate(col.id=1:nrow(.))) %>%
  dplyr::select(-geometry) %>%
  left_join(supp_2 %>%
              mutate(row.id=1:nrow(.))) %>%
  dplyr::select(-geometry)

county_intersection <- pi_2[, c(8, 13)]


supp_table <- supp %>%
  left_join(env_areas, by = "LOCALITY_ID") %>%
  left_join(county_intersection, by = "LOCALITY_ID") %>%
  rename("County" = "NAME")

write_csv(supp_table, "Data/supp_data/Supplementary_Table1.csv")
