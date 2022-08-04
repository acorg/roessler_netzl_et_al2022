rm(list = ls())
library(Racmacs)
library(ablandscapes)
library(tidyverse)
library(meantiter)
library(r3js)
library(htmlwidgets)
library(webshot2)
library(png)
library(grid)
library(gridExtra)
library(ggplot2)
library(patchwork)
set.seed(100)

source("./functions/landscape_functions.R")
source("./functions/remove_reactivity_bias.R")

landscape_dir <- file.path("figures", "landscapes")
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

single_exposure_sr_groups <- c("delta conv.", "alpha/alpha+E484K conv.","beta conv.","mRNA1273/mRNA1273","AZ/AZ","AZ/BNT","BNT/BNT","BA.1 conv." ,"BA.2 conv.", "WT conv.")

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
gmts_single <- plot_gmt_lndscps_from_map(map, sr_groups = c("WT conv.", "alpha/alpha+E484K conv.", "beta conv.", "delta conv.",
                                             "BA.1 conv.", "BA.2 conv."), adjust_reactivity_bias = F, remove_buttons = T)

gmts_single_plot <- plot_single_landscape_panel(gmts_single, label = "", label_x_pos = 1, label_y_pos = 0)
ggsave(file.path(landscape_dir, "gmts_single_exp.png"), gmts_single_plot, width = 4, height = 4, dpi = 300)

# plot double exposures
gmts_double_close <- plot_sr_group_lndscp_gmt(map, landscape_parameters,  sr_groups = c("mRNA1273/mRNA1273", "BNT/BNT", "AZ/BNT","AZ/AZ") , remove_buttons = T, adjust_reactivity_bias = F)

gmts_double_close_plot <- plot_single_landscape_panel(gmts_double_close, label = "", label_x_pos = 1, label_y_pos = 0)
ggsave(file.path(landscape_dir, "gmts_double_close_exp.png"), gmts_double_close_plot, width = 4, height = 4, dpi = 300)

# plot double exposures
gmts_double_dist <- plot_sr_group_lndscp_gmt(map, landscape_parameters,  sr_groups = c("BA.1 reinf.", "BA.2 reinf.") , remove_buttons = T, adjust_reactivity_bias = F)

gmts_double_dist_plot <- plot_single_landscape_panel(gmts_double_dist, label = "", label_x_pos = 1, label_y_pos = 0)
ggsave(file.path(landscape_dir, "gmts_double_dist_exp.png"), gmts_double_dist_plot, width = 4, height = 4, dpi = 300)

# plot triple exposures
gmts_triple_close <- plot_sr_group_lndscp_gmt(map, landscape_parameters, sr_groups =c("BNT/BNT/BNT","AZ/AZ+delta", "BNT/BNT+delta"),
  remove_buttons = T, adjust_reactivity_bias = T)

gmts_triple_close_plot <- plot_single_landscape_panel(gmts_triple_close, label = "", label_x_pos = 1, label_y_pos = 0)
ggsave(file.path(landscape_dir, "gmts_triple_close_exp.png"), gmts_triple_close_plot, width = 4, height = 4, dpi = 300)


# plot triple exposures
gmts_triple <- plot_sr_group_lndscp_gmt(map, landscape_parameters, sr_groups =c("Vacc+BA.1 reinf.",
                                                                 "Vacc+BA.2", "Vacc+BA.1") , remove_buttons = T, adjust_reactivity_bias = T)

gmts_triple_plot <- plot_single_landscape_panel(gmts_triple, label = "", label_x_pos = 1, label_y_pos = 0)
ggsave(file.path(landscape_dir, "gmts_triple_dist_exp.png"), gmts_triple_plot, width = 4, height = 4, dpi = 300)


# ------------------- Individual Landscape plots
labels <- LETTERS
names(labels) <- c("WT conv.", "alpha/alpha+E484K conv.", "beta conv.","delta conv.", "BA.1 conv.", "BA.2 conv.",
                                            "AZ/AZ", 'AZ/BNT', 'BNT/BNT', "mRNA1273/mRNA1273",
                                            "BA.1 reinf.", "BA.2 reinf.",
                                            "BNT/BNT/BNT","AZ/AZ+delta", "BNT/BNT+delta",
                                            "Vacc+BA.1","Vacc+BA.2", "Vacc+BA.1 reinf.")

plot_list <- list()
plot_list_no_letter <- list()

for(sg in c("WT conv.", "alpha/alpha+E484K conv.", "beta conv.","delta conv.", "BA.1 conv.", "BA.2 conv.",
                                            "AZ/AZ", 'AZ/BNT', 'BNT/BNT', "mRNA1273/mRNA1273")) {
  
  sg_label <- sg
  sg_label <- gsub("WT", "1st wave", sg_label)
  sg_label <- gsub("\\/alpha\\+E484K", "", sg_label)
  if(sg == "mRNA1273/mRNA1273"){
    sg_label <- "2xmRNA-1273"
  }
  
  if(sg_label %in% c("BA.1 conv.", "BA.2 conv.", "BA.5 conv.")) {
    sg_label <- gsub("conv.", "omicron conv.", sg_label)
  } 
  idvl_scape <- plot_sr_group_lndscp_from_map(map, sg, remove_buttons = TRUE, adjust_reactivity_bias = T)
  idvl_temp <- plot_single_landscape_panel(idvl_scape, label = labels[sg], label_x_pos = 0, label_y_pos = 9.5,
              label_size = 18,
              sr_group_label = sg_label, sr_group_size = 10, sr_group_y_pos = 0.5,
              show_border = TRUE)
  plot_list <- c(plot_list, list(idvl_temp))

  idvl_temp <- plot_single_landscape_panel(idvl_scape, label = "", label_x_pos = 0, label_y_pos = 9.5,
              label_size = 18,
              sr_group_label = sg_label, sr_group_size = 10, sr_group_y_pos = 0.5,
              show_border = TRUE)
  plot_list_no_letter <- c(plot_list_no_letter, list(idvl_temp))
  
}


for(sg in c("BA.1 reinf.", "BA.2 reinf.", "BNT/BNT/BNT","AZ/AZ+delta", "BNT/BNT+delta",
                                            "Vacc+BA.1","Vacc+BA.2", "Vacc+BA.1 reinf.")) {
 
  sg_label <- sg
  sg_label <- gsub("WT", "First wave", sg_label)

  if(dim(landscape_parameters[[sg]]$sr_cone_coords)[1]>0) {
    idvl_scape <- plot_sr_group_lndscp(map, landscape_parameters[[sg]]$log_titers, landscape_parameters[[sg]]$sr_cone_coords, landscape_parameters[[sg]]$colbases, landscape_parameters[[sg]]$slope, sg,
                              remove_buttons = TRUE, adjust_reactivity_bias =T)
    idvl_temp <- plot_single_landscape_panel(idvl_scape, label = labels[sg], label_x_pos = 0, label_y_pos = 9.5,
              label_size = 18,
              sr_group_label = sg_label, sr_group_size = 10, sr_group_y_pos = 0.5,
              show_border = TRUE)
  plot_list <- c(plot_list, list(idvl_temp))

  idvl_temp <- plot_single_landscape_panel(idvl_scape, label = "", label_x_pos = 0, label_y_pos = 9.5,
              label_size = 18,
              sr_group_label = sg_label, sr_group_size = 10, sr_group_y_pos = 0.5,
              show_border = TRUE)
  plot_list_no_letter <- c(plot_list_no_letter, list(idvl_temp))

  }
}

wrap_plots(plot_list, ncol = 4) -> idvl_scapes

ggsave(file.path(landscape_dir, "individual_landscapes_panel.png"), idvl_scapes, width = 16, height = 18, dpi = 300)

wrap_plots(plot_list_no_letter, ncol = 4) -> idvl_scapes

ggsave(file.path(landscape_dir, "individual_landscapes_panel_no_letter.png"), idvl_scapes, width = 16, height = 20, dpi = 300)

