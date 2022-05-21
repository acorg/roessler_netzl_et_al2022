# load required packages
rm(list = ls())
library(Racmacs)
library(stringr)

# set seed
set.seed(100)

# ------------------------------------------ FULL MAP --------------------------------------
full_map <- read.acmap("./data/maps/map-OmicronI+II+III-thresholded-full.ace")

full_map_p1_adj <- optimizeAgReactivity(full_map, fixed_ag_reactivities = c(0, 0, 0, -1, 0, 0, 0, 0), reoptimize = F)
full_map_p1_adj <- optimizeMap(full_map_p1_adj, number_of_dimensions = 2, number_of_optimizations = 1000, options = list(ignore_disconnected = TRUE))
full_map_p1_adj <- realignMap(full_map_p1_adj, full_map)

save.acmap(full_map_p1_adj, "./data/maps/map-OmicronI+II+III-thresholded-full-P1m1.ace")

##------------------------------- SINGLE EXPOSURE MAP ----------------------------------------
map_single <- read.acmap("./data/maps/map-OmicronI+II+III-thresholded-single_exposure.ace")

# optimize reactivity adjustment
single_map_p1_adj <- optimizeAgReactivity(map_single, fixed_ag_reactivities = c(0, 0, 0, -1, 0, 0, 0, 0), reoptimize = F)
single_map_p1_adj <- optimizeMap(single_map_p1_adj, number_of_dimensions = 2, number_of_optimizations = 1000, options = list(ignore_disconnected = TRUE))
single_map_p1_adj <- realignMap(single_map_p1_adj, map_single)

save.acmap(single_map_p1_adj, "./data/maps/map-OmicronI+II+III-thresholded-single_exposure-P1m1.ace")

