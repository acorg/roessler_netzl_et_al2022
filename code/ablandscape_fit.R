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

# function to fit slope and serum coordinates for multi exposure landscapes
fit_slope_all_sera <- function(
  pars,
  sr_coords,
  ag_coords,
  log_titers,
  colbase
) {
  
  if(length(colbase)>1) {
    pred_titers_all <- unlist(lapply(1:length(colbase), function(ser) {
      coords <- c(sr_coords[x,1], sr_coords[x,2])
      ag_distances <- as.matrix(dist(rbind(coords, ag_coords)))[1, -1]
      predicted_logtiters <- colbase[ser] - ag_distances*pars["slope"]
      return(predicted_logtiters)
    }))
  } else {
    coords <- c(sr_coords[1], sr_coords[2])
    ag_distances <- as.matrix(dist(rbind(coords, ag_coords)))[1, -1]
    pred_titers_all <- colbase - ag_distances*pars["slope"]
  }
  sum((log_titers - pred_titers_all)^2, na.rm = T)
  
}


# Read the base map
map <- read.acmap("./data/maps/map-OmicronI+II+III-thresholded-single_exposure-P1m1.ace")

# read the full map
map_orig <- read.acmap(paste0("./data/maps/map-OmicronI+II+III-thresholded-full-P1m1.ace"))

single_exposure_sr_groups <- c("delta conv.", "alpha/alpha+E484K conv.","beta conv.","mRNA1273/mRNA1273","AZ/AZ","AZ/BNT","BNT/BNT","BA.1 conv." ,"BA.2 conv.","BA.5 conv.", "WT conv.")

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

saveRDS(fit_results, "./data/landscape_fit/map-OmicronI+II+III-thresholded-full-P1m1-bias-no_1xVax+BA1.rds")
saveRDS(fit_results_no_bias, "./data/landscape_fit/map-OmicronI+II+III-thresholded-full-P1m1-bias_removed-no_1xVax+BA1.rds")

