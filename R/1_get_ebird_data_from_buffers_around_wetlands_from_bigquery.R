# get data from eBird
## packages
library(readr)
library(bigrquery)
library(dbplyr)
library(dplyr)
library(tidyr)
library(lubridate)
library(sf)
library(DBI)

# turn off spherical geometry!
# has a big speed increase when this is off. Like massive difference
# although I didn't time it out.
# this may be slightly less 'accurate', but some tests showed no difference
# some more info here: https://r-spatial.github.io/sf/articles/sf7.html
sf_use_s2(FALSE)

# read in study sites
study_sites <- st_read("Data/combined_shp_geojson/combined_polygons.geojson")


# create connection with online database
con <- DBI::dbConnect(bigrquery::bigquery(),
                      dataset= "ebird",
                      project="ebird-database",
                      billing="ebird-database")

# create ebird table
ebird <- tbl(con, 'ebird_qa_april_2022')

## extract data
# for all LOCALITY_ID from the combined dataset
localities <- ebird %>%
  dplyr::filter(OBSERVATION_DATE >= "2010-01-01") %>%
  dplyr::select(LOCALITY_ID, LATITUDE, LONGITUDE) %>%
  distinct() %>%
  collect(n=Inf)

# turn localities to sf
ebird_locations <- localities %>%
  st_as_sf(coords=c("LONGITUDE", "LATITUDE"), crs=4326)

# a function to filter obs within each polygon
poly_filter_obs_function <- function(poly){
  
  message(paste0("filtering observations for ", poly))
  
  locality <- study_sites %>%
    dplyr::filter(LOCALITY_ID==poly) %>%
    st_transform(crs=4326)
  
  lonlat2UTM = function(lonlat) {
    utm = (floor((lonlat[1] + 180) / 6) %% 60) + 1
    if(lonlat[2] > 0) {
      utm + 32600
    } else{
      utm + 32700
    }
  }
  
  coords <- st_coordinates(st_centroid(locality))
  EPSG_2_UTM <- lonlat2UTM(coords)
  # To see the UTM
  st_crs(EPSG_2_UTM)$proj4string
  
  locality_proj <- st_transform(st_as_sf(locality), EPSG_2_UTM)
  
  buffer_width_km <- 25
  
  buffered_locality <- locality_proj %>%
    st_buffer(dist = buffer_width_km * 1000)
  
  buffer_final <- buffered_locality %>%
    st_transform(crs=4326)

  # now filter the localities of all eBird locations
  locality_filter <- buffer_final %>%
    st_intersects(ebird_locations)
  
  # locations_within_buffer 
  locations_within_buffer <- localities %>%
    mutate(col.id=1:nrow(.)) %>%
    left_join(locality_filter %>%
                as.data.frame()) %>%
    dplyr::filter(row.id==1) %>%
    dplyr::select(-row.id, -col.id) %>%
    mutate(LOCALITY_ID_polygon=paste0(poly))
  
  return(locations_within_buffer)
  
}

all_localities_within_buffers <- bind_rows(lapply(unique(study_sites$LOCALITY_ID), poly_filter_obs_function))

length(unique(all_localities_within_buffers$LOCALITY_ID))

# now get the data using 75,000 localities at a time...
# because there is a character limit for a SQL query length


# now get the data....
dat_1 <- ebird %>%
  dplyr::filter(OBSERVATION_DATE >= "2010-01-01") %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, CATEGORY, BCR_CODE,
                TIME_OBSERVATIONS_STARTED, OBSERVER_ID, PROTOCOL_TYPE, DURATION_MINUTES, 
                EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE) %>%
  dplyr::filter(LOCALITY_ID %in% local(unique(all_localities_within_buffers$LOCALITY_ID)[1:75000])) %>%
  collect(n=Inf)

saveRDS(dat_1, "Data/eBird_buffer_data_chunks/chunk_1.RDS")

dat_2 <- ebird %>%
  dplyr::filter(OBSERVATION_DATE >= "2010-01-01") %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, CATEGORY, BCR_CODE,
                TIME_OBSERVATIONS_STARTED, OBSERVER_ID, PROTOCOL_TYPE, DURATION_MINUTES, 
                EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE) %>%
  dplyr::filter(LOCALITY_ID %in% local(unique(all_localities_within_buffers$LOCALITY_ID)[75001:150000])) %>%
  collect(n=Inf)

saveRDS(dat_2, "Data/eBird_buffer_data_chunks/chunk_2.RDS")

dat_3 <- ebird %>%
  dplyr::filter(OBSERVATION_DATE >= "2010-01-01") %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, CATEGORY, BCR_CODE,
                TIME_OBSERVATIONS_STARTED, OBSERVER_ID, PROTOCOL_TYPE, DURATION_MINUTES, 
                EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE) %>%
  dplyr::filter(LOCALITY_ID %in% local(unique(all_localities_within_buffers$LOCALITY_ID)[150001:225000])) %>%
  collect(n=Inf)

saveRDS(dat_3, "Data/eBird_buffer_data_chunks/chunk_3.RDS")

dat_4 <- ebird %>%
  dplyr::filter(OBSERVATION_DATE >= "2010-01-01") %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, CATEGORY, BCR_CODE,
                TIME_OBSERVATIONS_STARTED, OBSERVER_ID, PROTOCOL_TYPE, DURATION_MINUTES, 
                EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE) %>%
  dplyr::filter(LOCALITY_ID %in% local(unique(all_localities_within_buffers$LOCALITY_ID)[225001:300000])) %>%
  collect(n=Inf)

saveRDS(dat_4, "Data/eBird_buffer_data_chunks/chunk_4.RDS")

dat_5 <- ebird %>%
  dplyr::filter(OBSERVATION_DATE >= "2010-01-01") %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, CATEGORY, BCR_CODE,
                TIME_OBSERVATIONS_STARTED, OBSERVER_ID, PROTOCOL_TYPE, DURATION_MINUTES, 
                EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE) %>%
  dplyr::filter(LOCALITY_ID %in% local(unique(all_localities_within_buffers$LOCALITY_ID)[300001:375000])) %>%
  collect(n=Inf)

saveRDS(dat_5, "Data/eBird_buffer_data_chunks/chunk_5.RDS")

dat_6 <- ebird %>%
  dplyr::filter(OBSERVATION_DATE >= "2010-01-01") %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, CATEGORY, BCR_CODE,
                TIME_OBSERVATIONS_STARTED, OBSERVER_ID, PROTOCOL_TYPE, DURATION_MINUTES, 
                EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE) %>%
  dplyr::filter(LOCALITY_ID %in% local(unique(all_localities_within_buffers$LOCALITY_ID)[375001:450000])) %>%
  collect(n=Inf)

saveRDS(dat_6, "Data/eBird_buffer_data_chunks/chunk_6.RDS")

dat_7 <- ebird %>%
  dplyr::filter(OBSERVATION_DATE >= "2010-01-01") %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, CATEGORY, BCR_CODE,
                TIME_OBSERVATIONS_STARTED, OBSERVER_ID, PROTOCOL_TYPE, DURATION_MINUTES, 
                EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE) %>%
  dplyr::filter(LOCALITY_ID %in% local(unique(all_localities_within_buffers$LOCALITY_ID)[450001:525000])) %>%
  collect(n=Inf)

saveRDS(dat_7, "Data/eBird_buffer_data_chunks/chunk_7.RDS")

dat_8 <- ebird %>%
  dplyr::filter(OBSERVATION_DATE >= "2010-01-01") %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, CATEGORY, BCR_CODE,
                TIME_OBSERVATIONS_STARTED, OBSERVER_ID, PROTOCOL_TYPE, DURATION_MINUTES, 
                EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE) %>%
  dplyr::filter(LOCALITY_ID %in% local(unique(all_localities_within_buffers$LOCALITY_ID)[525001:600000])) %>%
  collect(n=Inf)

saveRDS(dat_8, "Data/eBird_buffer_data_chunks/chunk_8.RDS")

dat_9 <- ebird %>%
  dplyr::filter(OBSERVATION_DATE >= "2010-01-01") %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, CATEGORY, BCR_CODE,
                TIME_OBSERVATIONS_STARTED, OBSERVER_ID, PROTOCOL_TYPE, DURATION_MINUTES, 
                EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE) %>%
  dplyr::filter(LOCALITY_ID %in% local(unique(all_localities_within_buffers$LOCALITY_ID)[600001:675000])) %>%
  collect(n=Inf)

saveRDS(dat_9, "Data/eBird_buffer_data_chunks/chunk_9.RDS")

dat_10 <- ebird %>%
  dplyr::filter(OBSERVATION_DATE >= "2010-01-01") %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, CATEGORY, BCR_CODE,
                TIME_OBSERVATIONS_STARTED, OBSERVER_ID, PROTOCOL_TYPE, DURATION_MINUTES, 
                EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE) %>%
  dplyr::filter(LOCALITY_ID %in% local(unique(all_localities_within_buffers$LOCALITY_ID)[675001:750000])) %>%
  collect(n=Inf)

saveRDS(dat_10, "Data/eBird_buffer_data_chunks/chunk_10.RDS")

dat_11 <- ebird %>%
  dplyr::filter(OBSERVATION_DATE >= "2010-01-01") %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, CATEGORY, BCR_CODE,
                TIME_OBSERVATIONS_STARTED, OBSERVER_ID, PROTOCOL_TYPE, DURATION_MINUTES, 
                EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE) %>%
  dplyr::filter(LOCALITY_ID %in% local(unique(all_localities_within_buffers$LOCALITY_ID)[750001:774150])) %>%
  collect(n=Inf)

saveRDS(dat_11, "Data/eBird_buffer_data_chunks/chunk_11.RDS")