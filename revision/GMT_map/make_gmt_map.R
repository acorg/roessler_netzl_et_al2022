rm(list = ls())
library(Racmacs)
library(stringr)

# set seed
set.seed(100)

# load functions
source("./functions/map_functions.R")

# read in base map
map_single <- read.acmap("./data/maps/map-OmicronI+II+III-thresholded-single_exposure-P1m1.ace")

srGroup_colors <- read.csv(file = "./data/metadata/sr_group_colors.csv", header = TRUE, stringsAsFactors = FALSE, sep = ";", row.names = "SerumGroup")

sr_groups <- c("WT conv.", "delta conv.", "alpha/alpha+E484K conv.","beta conv.","mRNA1273/mRNA1273","AZ/AZ","AZ/BNT","BNT/BNT","BA.1 conv." ,"BA.2 conv.")

# do GMT map with <16 to 16/2
# calculate GMTs
log_table <- adjustedLogTiterTable(map_single)
log_table[round(2^log_table*10)==16] <- log2(8/10)
gmt_table <- sapply(sr_groups, function(x) {
    sub_table <- log_table[,as.character(srGroups(map_single)) == x]
    if(ncol(sub_table) > 1){
        rowMeans(sub_table, na.rm = T)
    } else {
        sub_table
    }
    
})

colnames(gmt_table) <- sr_groups

titer_gmt_table <- 2^gmt_table*10

gmt_map <- makeMap(titer_gmt_table, baseMap = map_single, nOptimisations = 1000, dilution_stepsize = 0, options = list(ignore_disconnected = TRUE))
srGroups(gmt_map) <- srNames(gmt_map)
srOutline(gmt_map) <- srGroup_colors[as.character(srGroups(gmt_map)),]
gmt_map <- apply_style(gmt_map)

gmt_map <- realignMap(gmt_map, map_single)
save.acmap(gmt_map, "revision/GMT_map/gmt_map_LOD2.ace")

