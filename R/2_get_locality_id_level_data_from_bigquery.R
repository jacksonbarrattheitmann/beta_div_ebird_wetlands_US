# get data from eBird
## This script will not run due to the API user specific key
## However, downloaded data are available in the 'Data' folder
## packages
library(readr)
library(bigrquery)
library(dbplyr)
library(dplyr)
library(tidyr)
library(lubridate)
library(sf)
library(DBI)

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
example_dat <- ebird %>%
  dplyr::filter(OBSERVATION_DATE >= "2010-01-01") %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, COMMON_NAME, SCIENTIFIC_NAME, OBSERVATION_COUNT,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, CATEGORY, BCR_CODE,
                TIME_OBSERVATIONS_STARTED, OBSERVER_ID, PROTOCOL_TYPE, DURATION_MINUTES, 
                EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER,
                COUNTRY, STATE_CODE) %>%
  dplyr::filter(LOCALITY_ID %in% local(study_sites$LOCALITY_ID)) %>%
  collect(n=Inf)

# filter to breeding season only
all_dat <- example_dat %>%
  mutate(MONTH=month(OBSERVATION_DATE, label=TRUE, abbr=TRUE))


# split data by month
# and write out each as an RDS
split_by_month_function <- function(month){
  
  temp <- all_dat %>% 
    dplyr::filter(MONTH==month)
  
  saveRDS(temp, paste0("Data/eBird_local_alpha_level/eBird_data_raw_", month, ".RDS"))
  
  
}

lapply(unique(all_dat$MONTH), split_by_month_function)
