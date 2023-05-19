# labelled map 
rm(list = ls())
library(Racmacs)

# read in map
full_map <- read.acmap("./data/maps/map-OmicronI+II+III-thresholded-single_exposure-P1m1.ace")
map <- read.acmap("revision/GMT_map/gmt_map_LOD2.ace")

xlim_zoom <- read.csv("./data/metadata/xlim_zoom.csv")$x
ylim_zoom <- read.csv("./data/metadata/ylim_zoom.csv")$x

xlim_no_zoom <- read.csv("./data/metadata/xlim_no_zoom.csv")$x
ylim_no_zoom <- read.csv("./data/metadata/ylim_no_zoom.csv")$x

# Setup plotting function
doplot <- function(map, xlims, ylims, text = TRUE) {
  
  # Setup the plot
  par(mar = rep(0.5, 4))
  
  agSize(map) <- 12
  # Plot the regular map
  srOutlineWidth(map) <- 1.2
  srOutline(map) <- adjustcolor(srOutline(map), alpha.f = 0.6)
  plot(map, xlim = xlims, 
       ylim =ylims, fill.alpha = 0.9)
  
  # Plot labels
  label_adjustments <- matrix(0, numAntigens(map), 2)
  rownames(label_adjustments) <- agNames(map)
  label_adjustments["B.1.351",] <- c(0.9, 0)
  label_adjustments["P.1.1",] <- c(-0.9, 0)
  label_adjustments["B.1.1.7+E484K",] <- c(0.9, -0.5)
  label_adjustments["BA.1",] <- c(-0.4, 0.7)
  label_adjustments["BA.2",] <- c(0, -0.6)
  label_adjustments["B.1.1.7",] <- c(0.3, -0.6)
  label_adjustments["D614G",] <- c(0, -0.5)
  label_adjustments["B.1.617.2",] <- c(0,-0.6)
  label_adjustments["BA.5",] <- c(0, 0.7)
  
  labels <- agNames(map)
  names(labels) <- agNames(map)
  labels["B.1.351"] <- "beta\n(B.1.351)"
  labels["P.1.1"] <- "gamma\n(P.1.1)"
  labels["BA.1"] <- "BA.1 omicron\n(B.1.1.529+BA.1)"
  labels["BA.2"] <- "BA.2 omicron\n(B.1.1.529+BA.2)"
  labels["B.1.617.2"] <- "delta\n(B.1.617.2)"
  labels["B.1.1.7"] <- "alpha\n(B.1.1.7)"
  labels["B.1.1.7+E484K"] <- "alpha + E484K\n(B.1.1.7+E484K)"
  labels["BA.5"] <- "BA.5 omicron\n(B.1.1.529+BA.5)"
  
  label_size <- rep(1, numAntigens(map))
  names(label_size) <- agNames(map)
  
  if(text){
         text(
    agCoords(map) + label_adjustments,
    cex = label_size,
    label = labels,
    font = 1
  )
  }
  
}

# make GMT map panel

png("revision/GMT_map/gmt_map_panel.png", width = 8, height = 6, units = 'in', res=300, pointsize = 18)
layout(matrix(c(1:4), ncol = 2, byrow = T))
par(oma=c(0, 0, 0, 0), mar=c(0.1, 0, 1, 0))

doplot(map, xlim_zoom, ylim_zoom, text = FALSE)
text(xlim_zoom[1]+0.4, ylim_zoom[2]-0.4, "A", cex = 1.2)

doplot(procrustesMap(map, full_map), xlim_zoom, ylim_zoom, text = FALSE)
text(xlim_zoom[1]+0.4, ylim_zoom[2]-0.4, "B", cex = 1.2)

plot(map, show_error = TRUE, xlim = xlim_zoom, ylim = ylim_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
       grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.3, plot_stress = TRUE)
text(xlim_zoom[1]+0.4, ylim_zoom[2]-0.4, "C", cex = 1.2)

plot(triangulationBlobs(relaxMap(map), stress_lim = 1, grid_spacing = 0.05), xlim = xlim_zoom, ylim = ylim_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
       grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.3, plot_stress = TRUE)
text(xlim_zoom[1]+0.4, ylim_zoom[2]-0.4, "D", cex = 1.2)

dev.off()
