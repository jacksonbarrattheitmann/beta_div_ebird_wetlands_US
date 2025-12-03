# Sampling polygons within the 25km buffer, that are the same area as 
# the GIW, then checking for checklist saturation equal to the GIW

library(readr)
library(tidyverse)
library(dplyr)
library(sf)

## Let's just try one to begin with
test <- readRDS("Data/ebird_data_by_ref_wetland_filt/L1009347.RDS")

test_stationary <- test %>%
  filter(PROTOCOL_TYPE == "Stationary")

num_check <- test_stationary %>%
  mutate(YEAR = str_sub(OBSERVATION_DATE, 1, 4)) %>%
  group_by(LOCALITY_ID, YEAR) %>%
  summarize(num_check = length(unique(SAMPLING_EVENT_IDENTIFIER)),
            ref_wetland = first(ref_wetland))

#the GIW
test_area <- st_read("Data/mapped_wetlands_geojson/L1009347.geojson")

# the 25km buffer around the GIW
test_area_25km <- st_buffer(test_area, dist = 25000)

# Now transform them both to UTM
test_area_utm <- transform_to_utm(test_area)

test_area_25km_utm <- st_transform(test_area_25km, crs=st_crs(test_area_utm))

random_points <- st_sample(test_area_25km_utm, size = 100, type = "random")


original_centroid <- st_centroid(st_union(test_area_utm))
orig_coords <- st_coordinates(original_centroid)

# For each random point, compute the translation offset
translated_polygons <- lapply(1:length(random_points), function(i) {
  target_point <- random_points[i]
  target_coords <- st_coordinates(target_point)
  
  # Offset
  dx <- target_coords[1] - orig_coords[1]
  dy <- target_coords[2] - orig_coords[2]
  
  # Apply translation
  moved_polygon <- st_geometry(test_area_utm) + c(dx, dy)
  st_sf(geometry = moved_polygon)
})

# Combine all into one sf object
random_polygons_sf <- do.call(rbind, translated_polygons)

# Need to set the crs
st_crs(random_polygons_sf) <- st_crs(test_area_utm)

#Now transform it back into the Long Lat projection to grab eBird data
random_polygons_wgs84 <- st_transform(random_polygons_sf, crs = 4326)

# Now we convert the eBird data in an sf object

test_sf <- test %>%
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"), crs = 4326)

# Check to make sure we have the same projection
st_crs(test_sf)
st_crs(random_polygons_wgs84)

# Looks good
# Now create a unique ID

random_polygons_wgs84 <- random_polygons_wgs84 %>%
  mutate(id = seq_len(nrow(random_polygons_wgs84)))


# Loop over polygons and filter checklists that intersect each
checklists_per_polygon <- map(
  random_polygons_wgs84$id,
  ~ {
    poly <- random_polygons_wgs84[random_polygons_wgs84$poly_id == .x, ]
    ebird_hits <- st_filter(test_sf, poly, .predicate = st_within)
    ebird_hits$polygon_id <- .x  # track which polygon it matched
    ebird_hits
  }
)

plot(st_geometry(random_polygons_wgs84), border = "red", col = NA, main = "eBird Points + Random Polygons")
plot(st_geometry(test_sf), pch = 20, col = "blue", add = TRUE)



# Read back in the RDS to make sure it looks good
test_polygons <- readRDS("Data/ebird_polygons_25km_buffer/ebird_checklists_by_polygon_test.rds") 
