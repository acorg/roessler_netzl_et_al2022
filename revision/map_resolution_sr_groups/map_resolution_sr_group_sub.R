library(Racmacs)
set.seed(100)
neut <- read.acmap('data/maps/map-OmicronI+II+III-thresholded-single_exposure-P1m1.ace')

xlim_no_zoom <- read.csv("./data/metadata/xlim_no_zoom.csv")$x
ylim_no_zoom <- read.csv("./data/metadata/ylim_no_zoom.csv")$x

neut_no_omi <- optimizeMap(removeSera(neut, srNames(neut)[srGroups(neut) %in% c('BA.1 conv.', 'BA.2 conv.')]), 2, 1000, options =  list(ignore_disconnected = TRUE))
neut_no_omi <- removeSera(neut_no_omi, c('C711_beta conv._B.1.351 conv._NA', 'C860_beta conv._B.1.351 conv._NA'))
neut_no_omi <- realignMap(neut_no_omi,neut)

map_positioned <- subsetMap(neut, sera = srNames(neut)[!(grepl("C711|C860|G347|G348|G353|G354|G355|G393", srNames(neut)))])

# select antigen in viewer and press "ctrl+c" to show connection lines
Racmacs::view(triangulationBlobs(neut_no_omi,grid_spacing = 0.05), options = list(xlim = xlim_no_zoom, ylim = ylim_no_zoom))
Racmacs::view(triangulationBlobs(relaxMap(map_positioned),grid_spacing = 0.05), options = list(xlim = xlim_no_zoom, ylim = ylim_no_zoom))
