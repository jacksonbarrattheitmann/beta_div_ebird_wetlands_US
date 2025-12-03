library(readr)
library(dplyr)
library(tidyverse)
library(data.table)
library(lattice)
library(geosphere)
library(purrr)
library(parallel)



# Use the ref list, to create 148 dfs with the gamma eBird data from each Locality ID in the buffer

# first read in the master ref list

combined_master_gamma_ref_wet <- readRDS(file = "Data/processed_buffer_chunks_with_ref_wetland/master_list_gamma_ref.RDS")


# creating the for loop to actually make 148 ref wetland data frames

# Might have to use cores overnight to run
##### registerDoParallel(detectCores())

process_wetland_data <- function(ref_wetland_id) {
  ind_wetland <- data.frame()
  
  for (chunk in 1:11) { 
    this_chunk <- readRDS(file = paste0("Data/eBird_buffer_data_chunks/chunk_",chunk,".RDS"))
    
    # should grab all the unique LOCALITY IDs associated with each ref wetland id
    gamma_hotspots_by_ref_wetland <- combined_master_gamma_ref_wet %>%
      filter(ref_wetland == ref_wetland_id) %>%
      select(LOCALITY_ID) %>%
      unique()
    
    # grabbing the actual eBird data from each chunk sequentially for the ref wetland
    eBird_data_by_ref <- this_chunk %>%
      filter(LOCALITY_ID %in% gamma_hotspots_by_ref_wetland$LOCALITY_ID)
    
    # binding each chunk output to a dataframe
    ind_wetland <- rbind(ind_wetland, eBird_data_by_ref)
  }
  
  # just adding a column with the ref wetland 
  ind_wetland <- ind_wetland %>%
    mutate(ref_wetland = ref_wetland_id)
  
  # saving an RDS for each unique reference wetland
  saveRDS(ind_wetland, file = paste0("Data/eBird_data_by_ref_wetland/", ref_wetland_id, ".RDS"))
}

# You can call the function like this:
unique_ref_wetlands <- unique(combined_master_gamma_ref_wet$ref_wetland)
for (wetland in unique_ref_wetlands) {
  process_wetland_data(wetland)
}

