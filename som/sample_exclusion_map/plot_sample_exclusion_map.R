rm(list = ls())
library(Racmacs)
set.seed(100)


original_map <- read.acmap("./data/maps/map-OmicronI+II+III-thresholded-single_exposure-P1m1.ace")

# subset map to exclude samples 
map <- subsetMap(original_map, sera = srNames(original_map)[!grepl("G776|G780", srNames(original_map))])
map <- optimizeMap(map, number_of_dimensions = 2, number_of_optimizations = 1000, options = list(ignore_disconnected = TRUE))
map <- realignMap(map, original_map)

save.acmap(map, "./data/maps/map-OmicronI+II+III-thresholded-single_exposure-P1m1-G776_G780_excluded.ace")

# read in subset map
map <- read.acmap("./data/maps/map-OmicronI+II+III-thresholded-single_exposure-P1m1-G776_G780_excluded.ace")

srFill(original_map)[grepl("G776|G780", srNames(original_map))] <- srOutline(original_map)[grepl("G776|G780", srNames(original_map))]

lims <- Racmacs:::mapPlotLims(original_map, sera = F)
lims$xlim <- read.csv("./data/metadata/xlim_zoom.csv")$x
lims$ylim <- read.csv("./data/metadata/ylim_zoom.csv")$x


do_proc_plot <- function(map, lims, letter) {
  
  # Setup the plot
  par(mar = c(0, 0.5, 0, 0.5))
  
  # Plot the regular map
  srOutlineWidth(map) <- 2
  srOutline(map) <- adjustcolor(srOutline(map), alpha.f = 0.3)
  plot(map, xlim = lims$xlim, ylim = lims$ylim, fill.alpha = 1)
  
  text(
    x = lims$xlim[1]+0.5, 
    y = lims$ylim[2]-0.6,
    cex = 2.5,
    label = letter,
    font = 1
  )
  
}


# Setup the plot
pdf(paste0("./som/sample_exclusion_map/ba2_cross_exclusion_map.pdf"), 8, 6)
layout(matrix(c(1:4), 2, 2, byrow = TRUE))

proc_map <- procrustesMap(map, original_map, optimization_number = 1, sera = FALSE)

do_proc_plot(original_map, lims, "A")
do_proc_plot(map, lims, "B")
do_proc_plot(proc_map, lims, "C")

dev.off()

png(paste0("./som/sample_exclusion_map/ba2_cross_exclusion_map.png"),  8, 6, units = "in", res = 300)
layout(matrix(c(1:4), 2, 2, byrow = TRUE))

proc_map <- procrustesMap(map, original_map, optimization_number = 1, sera = FALSE)

do_proc_plot(original_map, lims, "A")
do_proc_plot(map, lims, "B")
do_proc_plot(proc_map, lims, "C")

dev.off()
