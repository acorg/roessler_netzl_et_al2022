rm(list = ls())
library(Racmacs)
library(ablandscapes)
library(tidyverse)
library(meantiter)
library(r3js)
set.seed(100)

source("./functions/landscape_functions.R")
source("./functions/remove_reactivity_bias.R")


# Read the base map
map <- read.acmap("./data/maps/map-OmicronI+II+III-thresholded-single_exposure-P1m1-G776_G780_excluded.ace")


# read the full map
map_orig <- read.acmap(paste0("./data/maps/map-OmicronI+II+III-thresholded-full-P1m1.ace"))

# split srGroups of multiple infections to see how similar they are and if different groupings would make sense
sr_levels <- levels(srGroups(map_orig))
sr_levels <- c(sr_levels, "BA.2 G776", "BA.2 G780")

srGroups(map) <- factor(as.character(srGroups(map)), levels = sr_levels)
srGroups(map_orig) <- factor(as.character(srGroups(map_orig)), levels = sr_levels)

# assign the new srGroups to the specific sera
srGroups(map_orig)[grepl("G776", srNames(map_orig))] <- "BA.2 G776"
srGroups(map_orig)[grepl("G780", srNames(map_orig))] <- "BA.2 G780"

# get the sr groups to fit only those of interest
sr_groups_of_interest <- c("BA.2 G776", "BA.2 G780", "BA.2 reinf.", "Vacc+BA.2")
multi_exposure_sera <- srNames(map_orig)[(as.character(srGroups(map_orig)) %in% sr_groups_of_interest)]

# subset the map to only multi exposure sera
map_vax <- subsetMap(map_orig, sera = multi_exposure_sera)

# remove all only "*" sera from map (not titrated)
not_titrated <- apply(titerTable(map_vax), 2, function(x) all(x == "*"))
map_vax <- subsetMap(map_vax, sera = srNames(map_vax)[!(not_titrated)])

# do it for reactivity adjusted and unadjusted
log_titers <- adjustedLogTiterTable(map_vax)
colbase <- colBases(map_vax)
ag_coords <- unname(agCoords(map))


sr_groups <- unique(as.character(srGroups(map_vax)))
names(sr_groups) <- sr_groups

sr_group_colors <- data.frame("Serum.group" = as.character(srGroups(map_orig)),
                              "Color" = srOutline(map_orig)) %>% unique()


# # load parameters
# with reactivity bias removed
landscape_parameters <- readRDS("./som/landscapes_other_sr_groups/map-OmicronI+II+III-thresholded-full-P1m1-bias_removed-ba2_sr_groups.rds")

#a dd ba.2 conv to landscape parameters
# plot breakthrough infections omicron
sr_coords_ba2<- srCoords(map)[grepl("BA.2 conv.", rownames(srCoords(map))),]
sr_coords_ba2 <- sr_coords_ba2[!(is.na(sr_coords_ba2[,1])),]
sr_names <- rownames(sr_coords_ba2)
landscape_parameters[["BA.2 conv."]]$slope <- 1
landscape_parameters[["BA.2 conv."]]$log_titers <- adjustedLogTiterTable(map)[,sr_names]
landscape_parameters[["BA.2 conv."]]$colbases <- colBases(map)[srNames(map) %in% sr_names]
landscape_parameters[["BA.2 conv."]]$sr_cone_coords <- sr_coords_ba2

plot_sr_group_lndscp_gmt(map, landscape_parameters, sr_groups =c(sr_groups_of_interest, "BA.2 conv.") , remove_buttons = T, adjust_reactivity_bias = F)


for(sg in c(sr_groups_of_interest, "BA.2 conv.")) {
  
  
  if(dim(landscape_parameters[[sg]]$sr_cone_coords)[1]>0) {
    show(plot_sr_group_lndscp(map, fit_parameters_w_bias[[sg]]$log_titers, fit_parameters_w_bias[[sg]]$sr_cone_coords, fit_parameters_w_bias[[sg]]$colbases, fit_parameters_w_bias[[sg]]$slope, sg,
                              remove_buttons = TRUE, adjust_reactivity_bias = FALSE))
  }
}


