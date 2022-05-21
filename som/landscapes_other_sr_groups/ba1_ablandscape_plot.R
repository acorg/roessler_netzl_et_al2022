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
map <- read.acmap("./data/maps/map-OmicronI+II+III-thresholded-single_exposure-P1m1.ace")

# read the full map
map_orig <- read.acmap(paste0("./data/maps/map-OmicronI+II+III-thresholded-full-P1m1.ace"))

# split srGroups of multiple infections to see how similar they are and if different groupings would make sense
sr_levels <- levels(srGroups(map_orig))
sr_levels <- c(sr_levels, "1xVacc+BA.1", "2xVacc+BA.1", "3xVacc+BA.1", "B.1.617.2+BA.1", "D614G+BA.1", "B.1.1.7+BA.1")

srGroups(map) <- factor(as.character(srGroups(map)), levels = sr_levels)
srGroups(map_orig) <- factor(as.character(srGroups(map_orig)), levels = sr_levels)

# assign the new srGroups to the specific sera
ba1_vacc_ser <- srNames(map_orig)[srGroups(map_orig) == "Vacc+BA.1"]
nr_vacc <- unlist(lapply(ba1_vacc_ser, function(x) length(strsplit(x, "/")[[1]])))
ba1_sr_groups <- paste0(nr_vacc, "xVacc+BA.1")

srGroups(map_orig)[srNames(map_orig) %in% ba1_vacc_ser] <- ba1_sr_groups

ba1_inf_ser <- srNames(map_orig)[srGroups(map_orig) == "BA.1 reinf."]
srGroups(map_orig)[srNames(map_orig) %in% ba1_inf_ser] <- unlist(lapply(ba1_inf_ser, function(x) strsplit(x, "_")[[1]][3]))


# get the single exp sera
single_exposure_sr_groups <- c("delta conv.", "alpha/alpha+E484K conv.","beta conv.","mRNA1273/mRNA1273","AZ/AZ","AZ/BNT","BNT/BNT","BA.1 conv." ,"BA.2 conv.","WT conv.")

# get the sr groups to fit only multi exposure
multi_exposure_sera <- srNames(map_orig)[(as.character(srGroups(map_orig)) %in% c("Vacc+BA.1 reinf.", "2xVacc+BA.1", "3xVacc+BA.1", "B.1.617.2+BA.1", "D614G+BA.1", "B.1.1.7+BA.1"))]

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
fit_parameters_w_bias <- readRDS("./som/landscapes_other_sr_groups/map-OmicronI+II+III-thresholded-full-P1m1-bias-ba1_sr_groups.rds")

# with reactivity bias removed
landscape_parameters <- readRDS("./som/landscapes_other_sr_groups/map-OmicronI+II+III-thresholded-full-P1m1-bias_removed-ba1_sr_groups.rds")

# plot breakthrough infections omicron
plot_sr_group_lndscp_gmt(map, landscape_parameters, sr_groups =c("1xVacc+BA.1", "2xVacc+BA.1", "3xVacc+BA.1") , remove_buttons = F, adjust_reactivity_bias = F)

# plot omicron reinfections
plot_sr_group_lndscp_gmt(map, landscape_parameters, sr_groups =c("2xVacc+BA.1", "3xVacc+BA.1", "B.1.617.2+BA.1", "D614G+BA.1", "B.1.1.7+BA.1") , remove_buttons = T, adjust_reactivity_bias = F)

plot_sr_group_lndscp_gmt(map, landscape_parameters, sr_groups =c("3xVacc+BA.1", "Vacc+BA.1 reinf.") , remove_buttons = T, adjust_reactivity_bias = F)

# multi exposure, no need to adjust reactivity bias because it was adjusted already
plot_sr_group_lndscp_gmt(map, landscape_parameters, names(sr_groups), remove_buttons = FALSE, adjust_reactivity_bias = F)

for(sg in c("2xVacc+BA.1", "3xVacc+BA.1", "B.1.617.2+BA.1", "D614G+BA.1", "B.1.1.7+BA.1", "Vacc+BA.1 reinf.")) {
  
  
  if(dim(landscape_parameters[[sg]]$sr_cone_coords)[1]>0) {
    show(plot_sr_group_lndscp(map, fit_parameters_w_bias[[sg]]$log_titers, fit_parameters_w_bias[[sg]]$sr_cone_coords, fit_parameters_w_bias[[sg]]$colbases, fit_parameters_w_bias[[sg]]$slope, sg,
                              remove_buttons = TRUE, adjust_reactivity_bias = F))
  }
}

