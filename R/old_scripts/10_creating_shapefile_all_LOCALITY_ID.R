# Making shapefiles of all LOCALITY IDs
# within the 25km buffer of wetlands

library(tidyverse)
library(dplyr)
library(sf)


test_wet <- readRDS("Data/eBird_data_by_ref_wetland/L1009347.RDS")

master_gamma <- readRDS("Data/processed_buffer_chunks_with_ref_wetland/master_list_gamma_ref.RDS")

# first let's just get rid of all the checklists
# that are greater than 1KM in effort_distance
# or took over 240 minutes in Duration
test_wet_filt_check <- test_wet %>%
  filter(EFFORT_DISTANCE_KM <= 1 &
           DURATION_MINUTES < 240)

# Now we have to deal with group identifiers
get_group_checklists <- test_wet_filt_check %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, GROUP_IDENTIFIER) %>%
  distinct() %>%
  replace_na(list(GROUP_IDENTIFIER="not shared")) %>%
  dplyr::filter(GROUP_IDENTIFIER != "not shared") %>%
  group_by(GROUP_IDENTIFIER) %>%
  sample_n(1)

s_not_shared_checklists <- test_wet_filt_check %>%
  tidyr::replace_na(list(GROUP_IDENTIFIER="not shared")) %>%
  dplyr::filter(GROUP_IDENTIFIER == "not shared")

s_group_dat_trim <- test_wet_filt_check %>%
  filter(!duplicated(GROUP_IDENTIFIER))

s_final <- s_not_shared_checklists %>%
  bind_rows(s_group_dat_trim)

# summarize this data frame to get
# total checklists submitted at each individual LOCALITY_ID
# then filter to remove those with <25 checklists

wet_filt <- s_final %>%
  select(LOCALITY_ID, SAMPLING_EVENT_IDENTIFIER, LATITUDE, LONGITUDE) %>%
  group_by(LOCALITY_ID) %>%
  summarize(num_check = n_distinct(SAMPLING_EVENT_IDENTIFIER),
            LATITUDE = first(LATITUDE),
            LONGITUDE = first(LONGITUDE))

wet_filt_25 <- wet_filt %>%
  filter(num_check >= 25) 

# now we make into an sf_object
# create a buffer distance to characterize the site
wet_sf <- st_as_sf(wet_filt_25, coords = c("LONGITUDE", "LATITUDE"), crs = 4326)

wet_sf_1km <- st_buffer(wet_sf, dist = 250)

plot(wet_sf_1km)

wet_sf_1km <- wet_sf_1km %>%
  select(LOCALITY_ID, geometry)

# Now we can save this as a shapefile
write_sf(wet_sf_1km, "Data/1km_ebird_buffers/L1009347.shp")




# Creating a shapefile of ALL possible EBIRD HOTSPOTS at GAMMA scale

gamma_sf <- st_as_sf(master_gamma, coords = c("LONGITUDE", "LATITUDE"), crs = 4326)

gamma_sf <- gamma_sf %>%
  distinct(LOCALITY_ID)

st_write(gamma_sf, "Data/1km_ebird_buffers/all_shps.shp", driver = "ESRI Shapefile")

saveRDS(gamma_sf, "Data/1km_ebird_buffers/all_gammas_sf_point.RDS")

