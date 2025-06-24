# beta_div_ebird_wetlands_US
Spatial beta-diversity investigation of isolated wetlands in CONUS.

## R scripts

1_get_ebird_data_from_big_query: WILL NOT RUN. This script grabs the eBird data from CTC's cloud database with eBird data. It is provided in the data folder
under "Data/ebird_local_alpha_level/". 

4_creating_dataframe_for_just_GIWs_filtered: this script does a lot of the preliminary filtering of the data, including dealing with group submitted checklists, 'X' in teh abundnace column, and 'complete checklists'. The output is "Data/ebird_local_alpha_level/ebird_alpha_wetlands_raw". 

5_creating_env_matrix_for_alpha_ebird_data: this script appends Google Earth Engine data, wetland area, and EPA ecoregion level data into an env matrix that is 207 rows long (same number of GIWs). 

6_creating_site_species_matrix_GIWs_only: this is a preliminary script that creates the a 'wet_comm' dataframe object, that is a 207 row site x species matrix, with sites listed as .rownames. This species and abundance is aggregated across 73 checklists for each GIW, the smallest common denominator for submitted checklists across all sites. 

## Data
1. ebird_local_alpha_level:
   a. ebird_local_alpha_level_raw.RDS: this file is the cleaned and collated version of the other 9 RDS files in this folder. It is the one primarily used in the analysis, and for aggregating beta diversity among JUST GIWs. 

## Intermediate_data

1. all_wet_check_summarized.RDS: this file has the eBird data summarized by CHECKLIST and LOCALITY_ID for all 207 wetlands. Each LOCALITY ID has 73 rows, each representing a single checklist. This is used for the within GIW beta diversity analysis.
2. wet_comm_all_summarized.RDS: this file has the eBird data summarized by LOCALITY ID for all 207 wetlands. Each row is a single LOCALITY_ID with the summarized count data from all 73 checklists. This is a site x species matrix. 
3. wet_comm_ecoregion_10_summarized.RDS: this file is the same as above, but only contains 187 GIWS, that fall into ecoregions with MORE than 10 GIWs. This is used in the beta diversity analysis comaprions between ecoregions
4. env_ecoregions_10_GIWs.RDS: this file contains the ENVIRONMENTAL data for the above site x species matrix. 
