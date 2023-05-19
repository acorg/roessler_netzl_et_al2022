library(Racmacs)
set.seed(100)
neut <- read.acmap('data/maps/map-OmicronI+II+III-thresholded-single_exposure-P1m1.ace')

xlim_no_zoom <- read.csv("./data/metadata/xlim_no_zoom.csv")$x
ylim_no_zoom <- read.csv("./data/metadata/ylim_no_zoom.csv")$x

# subsetMap to only positioned sera
map_positioned <- subsetMap(neut, sera = srNames(neut)[!(grepl("C711|C860|G347|G348|G353|G354|G355|G393", srNames(neut)))])


png("revision/error_triangulation_blobs/error_triangulation.png", width = 6.5, height = 3, units = 'in', res=300, pointsize = 18)
layout(matrix(c(1:2), ncol = 2, byrow = T))
par(oma=c(0, 0, 0, 0), mar=c(0.1, 0, 1, 0))
plot(map_positioned, show_error = TRUE, xlim = xlim_no_zoom, ylim = ylim_no_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
       grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.3, plot_stress = TRUE)
text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, "A", cex = 1.2)

plot(triangulationBlobs(relaxMap(map_positioned), stress_lim = 1, grid_spacing = 0.05), xlim = xlim_no_zoom, ylim = ylim_no_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
       grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.3, plot_stress = TRUE)
text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, "B", cex = 1.2)

dev.off()
