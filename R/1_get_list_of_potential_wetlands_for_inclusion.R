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

# create connection with online database
con <- DBI::dbConnect(bigrquery::bigquery(),
                      dataset= "ebird",
                      project="ebird-database",
                      billing="ebird-database")

# create ebird table
ebird <- tbl(con, 'ebird_qa_april_2022')

## extract data
# with strict filters
## extract data
example_dat <- ebird %>%
  dplyr::filter(OBSERVATION_DATE >= "2015-01-01") %>%
  dplyr::select(SAMPLING_EVENT_IDENTIFIER, LOCALITY, LOCALITY_TYPE, COUNTRY,
                LOCALITY_ID, LATITUDE, LONGITUDE, OBSERVATION_DATE, BCR_CODE,
                TIME_OBSERVATIONS_STARTED, OBSERVER_ID, PROTOCOL_TYPE, DURATION_MINUTES, 
                EFFORT_DISTANCE_KM, EFFORT_AREA_HA, NUMBER_OBSERVERS, GROUP_IDENTIFIER, STATE_CODE) %>%
  dplyr::filter(COUNTRY=="United States") %>%
  group_by(LOCALITY, LOCALITY_TYPE, LOCALITY_ID, LATITUDE, LONGITUDE, STATE_CODE) %>%
  summarize(Number_of_checklists=n_distinct(SAMPLING_EVENT_IDENTIFIER)) %>%
  dplyr::filter(Number_of_checklists>=10) %>%
  collect(n=Inf)


# now look to see how many hotspots have the word wetland or wetlands in it
wetlands <- example_dat %>%
  mutate(has_wetland=stringr::str_detect(LOCALITY, "wetland")) %>%
  mutate(has_wetland2=stringr::str_detect(LOCALITY, "Wetland")) %>%
  dplyr::filter(has_wetland=="TRUE" | has_wetland2=="TRUE") %>%
  arrange(desc(Number_of_checklists)) %>%
  dplyr::select(-has_wetland, -has_wetland2) %>%
  dplyr::filter(! STATE_CODE %in% c("US-AK", "US_HI"))

write_csv(wetlands, "Data/potential_study_sites_list/potential_wetland_ebird_locations.csv")




