# labelled map 
rm(list = ls())
library(Racmacs)

# read in map
map <- read.acmap("./data/maps/map-OmicronI+II+III-thresholded-single_exposure-P1m1.ace")

xlim_zoom <- read.csv("./data/metadata/xlim_zoom.csv")$x
ylim_zoom <- read.csv("./data/metadata/ylim_zoom.csv")$x

xlim_no_zoom <- read.csv("./data/metadata/xlim_no_zoom.csv")$x
ylim_no_zoom <- read.csv("./data/metadata/ylim_no_zoom.csv")$x

# Setup plotting function
doplot <- function(map, xlims, ylims) {
  
  # Setup the plot
  par(mar = rep(0.5, 4))
  
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
  
  text(
    agCoords(map) + label_adjustments,
    cex = label_size,
    label = labels,
    font = 1
  )
  
}

pdf("figures/labelled_map/labelled_map_no_zoom.pdf", 7, 7)
par(mar = rep(0.5, 4))
doplot(map, xlim_no_zoom, ylim_no_zoom)
dev.off()

png("figures/labelled_map/labelled_map_no_zoom.png", 7, 7, units = 'in', res=300, pointsize = 12)
par(mar = rep(0.5, 4))
doplot(map, xlim_no_zoom, ylim_no_zoom)
dev.off()


pdf("figures/labelled_map/labelled_map_zoom.pdf", 5, 4)
par(mar = rep(0.5, 4))
doplot(map, xlim_zoom, ylim_zoom)
dev.off()

png("figures/labelled_map/labelled_map_zoom.png", 5, 4, units = 'in', res=300, pointsize = 12)
par(mar = rep(0.5, 4))
doplot(map, xlim_zoom, ylim_zoom)
dev.off()
