library(sf)
library(dplyr)
library(readr)
library(tidyverse)
library(terra)

gdb_path <- "Data/global_NFWs/Global_NFWs_NorthAmerica.gdb/"

st_layers(gdb_path)

layer_name <- "Global_NFWs_7040000010"

sf_layer <- read_sf(gdb_path, layer = layer_name)

raster_layers <- terra::rast(gdb_path)
