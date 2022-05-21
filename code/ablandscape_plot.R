rm(list = ls())
library(Racmacs)
library(ablandscapes)
library(tidyverse)
library(meantiter)
library(r3js)
set.seed(100)

source("./functions/landscape_functions.R")
source("./functions/remove_reactivity_bias.R")


# function to add slope 1 landscapes to landscape parameters
add_sr_group_parameters <- function(landscape_parameters, sr_groups, map){
  
  for(sr_group in sr_groups){
    sr_coords<- srCoords(map)[grepl(sr_group, rownames(srCoords(map))),]
    sr_coords <- sr_coords[!(is.na(sr_coords[,1])),]
    sr_names <- rownames(sr_coords)
    landscape_parameters[[sr_group]]$slope <- 1
    landscape_parameters[[sr_group]]$log_titers <- adjustedLogTiterTable(map)[,sr_names]
    landscape_parameters[[sr_group]]$colbases <- colBases(map)[srNames(map) %in% sr_names]
    landscape_parameters[[sr_group]]$sr_cone_coords <- sr_coords
  }
  
  
  return(landscape_parameters)
  
}

# Read the base map
map <- read.acmap("./data/maps/map-OmicronI+II+III-thresholded-single_exposure-P1m1.ace")

# read the full map
map_orig <- read.acmap(paste0("./data/maps/map-OmicronI+II+III-thresholded-full-P1m1.ace"))

single_exposure_sr_groups <- c("delta conv.", "alpha/alpha+E484K conv.","beta conv.","mRNA1273/mRNA1273","AZ/AZ","AZ/BNT","BNT/BNT","BA.1 conv." ,"BA.2 conv.","WT conv.")

single_exposure_sr <- srNames(map_orig)[as.character(srGroups(map_orig)) %in% single_exposure_sr_groups]

# get the sr groups to fit only multi exposure
multi_exposure_sera <- srNames(map_orig)[!(as.character(srGroups(map_orig)) %in% single_exposure_sr_groups)]

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
# reactivity bias not removed
fit_parameters_w_bias <- readRDS("./data/landscape_fit/map-OmicronI+II+III-thresholded-full-P1m1-bias-no_1xVax+BA1.rds")

# with reactivity bias removed
landscape_parameters <- readRDS("./data/landscape_fit/map-OmicronI+II+III-thresholded-full-P1m1-bias_removed-no_1xVax+BA1.rds")

# add 2x vax landscape parameters
landscape_parameters <- add_sr_group_parameters(landscape_parameters, sr_groups = c("mRNA1273/mRNA1273", "BNT/BNT", "AZ/BNT","AZ/AZ"), map = map)

# plot unvaccinated convalescent landscapes
plot_gmt_lndscps_from_map(map, sr_groups = c("WT conv.", "alpha/alpha+E484K conv.", "delta conv.", "beta conv.",
                                             "BA.1 conv.", "BA.2 conv."), adjust_reactivity_bias = F, remove_buttons = T)

# plot unvaccinated convalescent landscapes
plot_gmt_lndscps_from_map(map, sr_groups = c("WT conv.", "alpha/alpha+E484K conv.", "delta conv.", "beta conv.",
                                             "BA.1 conv.", "BA.2 conv."), adjust_reactivity_bias = T, remove_buttons = T)


# plot double exposures
plot_sr_group_lndscp_gmt(map, landscape_parameters,  sr_groups = c("mRNA1273/mRNA1273", "BNT/BNT", "AZ/BNT","AZ/AZ", "BA.1 reinf.", "BA.2 reinf.") , remove_buttons = T, adjust_reactivity_bias = T)
# plot double exposures
plot_sr_group_lndscp_gmt(map, landscape_parameters,  sr_groups = c("mRNA1273/mRNA1273", "BNT/BNT", "AZ/BNT","AZ/AZ", "BA.1 reinf.", "BA.2 reinf.") , remove_buttons = F, adjust_reactivity_bias = F)


# plot triple exposures
plot_sr_group_lndscp_gmt(map, landscape_parameters, sr_groups =c("AZ/AZ+delta", "BNT/BNT+delta","Vacc+BA.1 reinf.",
                                                                 "Vacc+BA.2", "Vacc+BA.1") , remove_buttons = T, adjust_reactivity_bias = T)


plot_gmt_lndscps_from_map(map, adjust_reactivity_bias = T, remove_buttons = T)
plot_gmt_lndscps_from_map(map, adjust_reactivity_bias = F)


for(sg in unique(as.character(srGroups(map)))) {
  print(sg)
  # ony adjusts GMT impuluses. does not change shape of landscapes
  show(plot_sr_group_lndscp_from_map(map, sg, remove_buttons = TRUE, adjust_reactivity_bias = T))
}

show(plot_sr_group_lndscp_from_map(map, "BA.1 conv.", remove_buttons = TRUE, adjust_reactivity_bias = T))
show(plot_sr_group_lndscp_from_map(map, "beta conv.", remove_buttons = F, adjust_reactivity_bias = F))

# multi exposure, no need to adjust reactivity bias because it was adjusted already
plot_sr_group_lndscp_gmt(map, landscape_parameters, names(sr_groups), remove_buttons = FALSE, adjust_reactivity_bias = F)

for(sg in names(landscape_parameters)) {
 
  print(sg)
  if(dim(landscape_parameters[[sg]]$sr_cone_coords)[1]>0) {
    show(plot_sr_group_lndscp(map, landscape_parameters[[sg]]$log_titers, landscape_parameters[[sg]]$sr_cone_coords, landscape_parameters[[sg]]$colbases, landscape_parameters[[sg]]$slope, sg,
                              remove_buttons = TRUE, adjust_reactivity_bias =T))
  }
  
  # if(dim(landscape_parameters[[sg]]$sr_cone_coords)[1]>0) {
  #   show(plot_sr_group_lndscp(map, fit_parameters_w_bias[[sg]]$log_titers, fit_parameters_w_bias[[sg]]$sr_cone_coords, fit_parameters_w_bias[[sg]]$colbases, fit_parameters_w_bias[[sg]]$slope, sg,
  #                             remove_buttons = TRUE, adjust_reactivity_bias = FALSE))
  # }
}

for(sg in c("Vacc+BA.1", "3xVacc+BA.2")) { #} names(sr_groups)) {
  
  print(sg)
  if(dim(landscape_parameters[[sg]]$sr_cone_coords)[1]>0) {
    show(plot_sr_group_lndscp(map, fit_parameters_w_bias[[sg]]$log_titers, fit_parameters_w_bias[[sg]]$sr_cone_coords, fit_parameters_w_bias[[sg]]$colbases, fit_parameters_w_bias[[sg]]$slope, sg,
                              remove_buttons = TRUE, adjust_reactivity_bias = FALSE))
  }
}

srNames(map_vax)[srGroups(map_vax) == "Vacc+BA.2"]

