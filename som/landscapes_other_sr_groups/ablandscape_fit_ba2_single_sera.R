#setup page and load metadata
rm(list = ls())
library(Racmacs)
library(tidyverse)
library(meantiter)
set.seed(100)

source("./functions/remove_reactivity_bias.R")

# function to fit slope and serum coordinates for multi exposure landscapes
fit_cone_all_sera <- function(
  pars,
  ag_coords,
  log_titers,
  colbase
) {
  
  if(length(colbase)>1) {
    pred_titers_all <- unlist(lapply(1:length(colbase), function(ser) {
      coords <- c(pars[paste0("x",ser)], pars[paste0("y",ser)])
      ag_distances <- as.matrix(dist(rbind(coords, ag_coords)))[1, -1]
      predicted_logtiters <- colbase[ser] - ag_distances*pars["slope"]
      return(predicted_logtiters)
    }))
  } else {
    coords <- c(pars["x"], pars["y"])
    ag_distances <- as.matrix(dist(rbind(coords, ag_coords)))[1, -1]
    pred_titers_all <- colbase - ag_distances*pars["slope"]
  }
  sum((log_titers - pred_titers_all)^2, na.rm = T)
  
}


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

# fit for each serum group level and remove reactivity bias
fit_results_no_bias <- lapply(as.character(unique(srGroups(map_vax))), function(x) {
  log <-as.matrix(log_titers[,as.character(srGroups(map_vax)) == x])
  
  unadj_log <- log
  unadj_colb <- colbase[as.character(srGroups(map_vax)) == x]
  # remove reactivity bias for each serum
  log <- remove_reactivity_bias_logtiter(log)
  colb <- apply(log, 2, function(x) max(x, na.rm = T))
  
  fit <-  optim(
    par = c(
      x = ag_coords[apply(log, 2, which.max),1],
      y = ag_coords[apply(log, 2, which.max),2],
      slope = 1
    ),
    fn = fit_cone_all_sera,
    method = "L-BFGS-B",
    ag_coords = ag_coords,
    log_titers = log,
    colbase = colb,
    control = list(maxit = 500)
  )
  
  pars <- fit[["par"]]
  fit$slope <- pars["slope"]
  
  fit$sr_cone_coords <- as.matrix(cbind(pars[grepl("x", names(pars))], pars[grepl("y", names(pars))]))
  
  fit$colbases_adj <- colb
  
  fit$log_titers_adj <- log
  
  fit$colbases <- unadj_colb
  
  fit$log_titers <- unadj_log
  
  return(fit)
})

# fit for each serum group level without removing reactivity bias
fit_results <- lapply(as.character(unique(srGroups(map_vax))), function(x) {
  log <-as.matrix(log_titers[,as.character(srGroups(map_vax)) == x])
  colb <- colbase[as.character(srGroups(map_vax)) == x]
  
  fit <-  optim(
    par = c(
      x = ag_coords[apply(log, 2, which.max),1],
      y = ag_coords[apply(log, 2, which.max),2],
      slope = 1
    ),
    fn = fit_cone_all_sera,
    method = "L-BFGS-B",
    ag_coords = ag_coords,
    log_titers = log,
    colbase = colb,
    control = list(maxit = 500)
  )
  
  pars <- fit[["par"]]
  fit$slope <- pars["slope"]
  
  fit$sr_cone_coords <- as.matrix(cbind(pars[grepl("x", names(pars))], pars[grepl("y", names(pars))]))
  
  fit$colbases <- colb
  
  fit$log_titers <- log
  
  return(fit)
})

names(fit_results_no_bias) <- as.character(unique(srGroups(map_vax)))
names(fit_results) <- as.character(unique(srGroups(map_vax)))

saveRDS(fit_results, "./som/landscapes_other_sr_groups/map-OmicronI+II+III-thresholded-full-P1m1-bias-ba2_sr_groups.rds")
saveRDS(fit_results_no_bias, "./som/landscapes_other_sr_groups/map-OmicronI+II+III-thresholded-full-P1m1-bias_removed-ba2_sr_groups.rds")
