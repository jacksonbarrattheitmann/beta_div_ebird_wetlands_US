# Map of GIWs w/ ecoregion

library(sf)
library(dplyr)
library(readr)
library(tidyverse)
library(ggplot2)
library(maps)
library(tigris)
options(tigris_use_cache = TRUE)
library(viridis)
# read in the ecoregion spatial data

eco <- read_sf("Data/earth_engine_env_data/epa_level1_eco/NA_CEC_Eco_Level1.shp")

# Get the full USA (state boundaries)
usa_states <- states(cb = TRUE)

# Filter to 50 states only if needed
usa_states <- usa_states %>%
  filter(!STUSPS %in% c("AS", "GU", "MP", "PR", "VI", "AK", "HI"))  # Remove territories and HI and AK

usa_outline <- usa_states %>% 
  st_union() %>% 
  st_as_sf()


# Ensure CRS matches
eco <- st_transform(eco, st_crs(usa_outline))

# Make geomertires valid
eco <- st_make_valid(eco) 

# Clip
ecoregions_us <- st_intersection(eco, usa_outline)

# Union (dissolve) by ecoregion column
ecoregions_dissolved <- ecoregions_us %>%
  group_by(NA_L1NAME) %>%
  summarise(geometry = st_union(geometry), .groups = "drop") %>%
  st_as_sf()

# Simplify geometry (tolerance is in map units, e.g., degrees if unprojected)
ecoregions_simplified <- st_simplify(ecoregions_dissolved, dTolerance = 0.10, preserveTopology = TRUE)

ecoregions_simplified <- ecoregions_simplified %>%
  filter(!NA_L1NAME == "WATER")

# Get unique ecoregions and assign colors
ecoregion_names <- unique(ecoregions_simplified$NA_L1NAME)
n_colors <- length(ecoregion_names)
colors <- viridis(n_colors, option = "C") 

# Match colors to ecoregion factor levels
ecoregion_factor <- factor(ecoregions_simplified$NA_L1NAME)
fill_colors <- colors[as.numeric(ecoregion_factor)]

# Now need to load in the GIW spatial data to overlay on the plot
giw_spat <- readRDS("Data/eBird_local_alpha_level/ebird_alpha_wetlands_raw.RDS") %>%
  select(LOCALITY_ID, LATITUDE, LONGITUDE) %>%
  distinct(LOCALITY_ID, .keep_all = TRUE)

wetlands_sf <- st_as_sf(giw_spat, 
                        coords = c("LONGITUDE", "LATITUDE"), 
                        crs = 4326)

wetlands_sf <- st_transform(wetlands_sf, st_crs(ecoregions_simplified))

# Spatial join: assign each wetland the ecoregion it's in
wetlands_joined <- st_join(wetlands_sf, ecoregions_simplified[, "NA_L1NAME"])


# Factor to match plotting order
wetland_ecoregion_factor <- factor(wetlands_joined$NA_L1NAME, 
                                   levels = levels(ecoregion_factor))

# Use same colors as in the map (from earlier)
point_border_colors <- colors[as.numeric(wetland_ecoregion_factor)]


# Plot the simplified ecoregions
plot(ecoregions_simplified["NA_L1NAME"], 
     main = "",
     col = fill_colors,       # Fill by ecoregion
     border = "black", 
     reset = FALSE, 
     key.pos = NULL)          # Disable default legend

plot(wetlands_joined$geometry,
     add = TRUE,
     pch = 21,               # circle with fill and border
     bg = point_border_colors,           # fill color by ecoregion
     col = "darkgray",  # border color
     lwd = 2,              # thicker outline
     cex = 1.5)


# Add a custom legend at the bottom
legend("bottom", 
       legend = levels(ecoregion_factor), 
       fill = colors,
       ncol = 3,              # Spread across columns
       bty = "n",             # No border
       cex = 0.75,             # Smaller text
       inset = c(0, -0.05), 
       xpd = TRUE)   




