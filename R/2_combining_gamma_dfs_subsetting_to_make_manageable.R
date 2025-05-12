library(readr)
library(dplyr)
library(tidyverse)
library(data.table)
library(lattice)
library(geosphere)
library(purrr)

# files <- list.files(path = "Data/eBird_buffer_data_chunks/", pattern = "\\.RDS", full.names = TRUE)
# big_r <- lapply(files, readRDS)
# big_s <- bind_rows(big_r, .id = "column_label")
# big_s <- as.data.frame(big_s)

# load in JUST wetlands
file_list <- list.files(path = "Data/eBird_local_alpha_level/", pattern = "\\.RDS", full.names = TRUE)
r <- lapply(file_list, readRDS)
s <- bind_rows(r, .id = "column_label")
s <- as.data.frame(s)

# list of the unique locality IDs that are wetlands, with the geometry (LAT, LONG)

wetlands_ll <- s %>%
  group_by(LOCALITY_ID) %>%
  summarise(LATITUDE = unique(LATITUDE), LONGITUDE = unique(LONGITUDE)) %>%
  select(LOCALITY_ID, LONGITUDE, LATITUDE)


for (chunk in 1:11){ 
  # Loading in the first gamma chunk to test
  # Let's just try loading in 1 chunk
  
  gamma_chunk <- readRDS(file = paste0("Data/eBird_buffer_data_chunks/chunk_",chunk,".RDS"))
  
  # Generate data frame with LOCALITY_ID, LAT, LONG, for the Gamma chunk
  
  gamma_ll <- gamma_chunk %>%
    group_by(LOCALITY_ID) %>%
    summarise(LATITUDE = unique(LATITUDE), LONGITUDE = unique(LONGITUDE)) %>%
    select(LOCALITY_ID, LONGITUDE, LATITUDE)

  gamma <- gamma_ll
  
  gamma_chunk_wetland_final <- data.frame()
# for loop going through rownumber in 1:length(wetlands_ll$LOCALITY_ID)
for(rownumber in 1:length(wetlands_ll$LOCALITY_ID)){
       
       test_dist_matrix <- as.data.frame(t(distm(wetlands_ll[rownumber, 2:3], gamma[ , 2:3]))) %>%
         mutate(LOCALITY_ID = gamma$LOCALITY_ID) %>%
         rename(dist = V1)
       
       # checking to see if gamma hotspots are within buffer to first wetland LOCALITY_ID
       
       within_buffer <- test_dist_matrix %>%
         filter(dist < 25000)
       
       
       # creating the LOCALITY_ID sample
       all_ids <- as.vector(unique(within_buffer$LOCALITY_ID))
       
       gamma_chunk_wetland <- gamma %>%
         filter(LOCALITY_ID %in% all_ids) %>%
         mutate(ref_wetland = wetlands_ll$LOCALITY_ID[rownumber])
       
       # this is not getting the minimum distance to closest wetland, just taking it 
       # based on the list of wetland_ids (should fix at some point)
       
       gamma_chunk_wetland_final <- rbind(gamma_chunk_wetland_final, gamma_chunk_wetland)
  }

  # Save the RDS
  saveRDS(gamma_chunk_wetland_final, file = paste0("Data/processed_buffer_chunks_with_ref_wetland/processed_gamma_chunk",chunk,".RDS"))
}

 
# combine these into 1 master file with all the ref wetlands

combined_master_gamma_ref_wet <- map_dfr(1:11, .f = function(chunk){
  readRDS(file = paste0("Data/processed_buffer_chunks_with_ref_wetland/processed_gamma_chunk",chunk,".RDS"))
}
  )

saveRDS(combined_master_gamma_ref_wet, file = "Data/processed_buffer_chunks_with_ref_wetland/master_list_gamma_ref.RDS")




