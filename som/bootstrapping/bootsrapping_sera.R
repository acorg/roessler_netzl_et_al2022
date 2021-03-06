library(Racmacs)
set.seed(100)
neut <- read.acmap('data/maps/map-OmicronI+II+III-thresholded-single_exposure-P1m1.ace')

xlim_no_zoom <- read.csv("./data/metadata/xlim_no_zoom.csv")$x
ylim_no_zoom <- read.csv("./data/metadata/ylim_no_zoom.csv")$x-1

labels <- data.frame(
  row.names = c('mRNA1273/mRNA1273', 'AZ/AZ', 'AZ/BNT', 'BNT/BNT',"WT conv.", 'alpha/alpha+E484K conv.', 'beta conv.', 'delta conv.', 'BA.1 conv.', 'BA.2 conv.'),
  val = c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J')
)


png("som/bootstrapping/bootstrapping-sera.png", width = 8, height = 8, units = 'in', res=300, pointsize = 18)
layout(matrix(c(1:12), ncol = 3, byrow = T))
par(oma=c(0, 0, 0, 0), mar=c(0.1, 0, 1, 0))

for(srGroup in c('mRNA1273/mRNA1273', 'AZ/AZ', 'AZ/BNT', 'BNT/BNT',"WT conv.", 'alpha/alpha+E484K conv.', 'beta conv.', 'delta conv.', 'BA.1 conv.', 'BA.2 conv.')){
    
    print(srGroup)
  
  newMap <- removeSera(neut, srNames(neut)[srGroups(neut) == srGroup])
  newMap <- optimizeMap(
    map                     = newMap,
    number_of_dimensions    = 2,
    number_of_optimizations = 500,
    minimum_column_basis    = "none",
    options =  list(ignore_disconnected = TRUE)
  )
  
  newMap <- realignMap(newMap, neut)
  
  save_text <- strsplit(srGroup, "/")[[1]][1]
  if(srGroup == "AZ/BNT"){
    save_text <- "AZ_BNT"
  }
  
  # uncomment below if you want to save map
  #save.acmap(map = newMap, filename = paste0("./som/bootstrapping/wo_",save_text,".ace"))
  srOutlineWidth(newMap) <- 1
  
  p <- procrustesMap(newMap, neut, sera = FALSE)
  
  title_text <- srGroup
  if(title_text == "mRNA1273/mRNA1273") {
    title_text <- "mRNA-1273/mRNA-1273"
  } else if (title_text == "WT conv.") {
    title_text <- "first wave conv."
  } else if (title_text == "BA.1 conv.") {
    title_text <- "BA.1 omicron conv."
  } else if (title_text == "BA.2 conv.") {
    title_text <- "BA.2 omicron conv."
  } else if (title_text == "alpha/alpha+E484K conv.") {
    title_text <- "alpha conv."
  }
  
  plot(p, xlim = xlim_no_zoom, ylim = ylim_no_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
       grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.3, plot_stress = TRUE)
  title(main = paste(title_text, 'sera', sep=' '), cex.main=0.7, line = 0.2)
  text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, labels[srGroup, ], cex = 1.2)
}

dev.off()

stop()

# run below if maps were saved above
# when maps already exist
png("som/bootstrapping/bootstrapping-sera.png", width = 10, height = 6, units = 'in', res=300, pointsize = 18)
layout(matrix(c(1:12), ncol = 4, byrow = T))
par(oma=c(0, 0, 0, 0), mar=c(0.1, 0, 1, 0))

for(srGroup in c('mRNA1273/mRNA1273', 'AZ/AZ', 'AZ/BNT', 'BNT/BNT',"WT conv.", 'alpha/alpha+E484K conv.', 'beta conv.', 'delta conv.', 'BA.1 conv.', 'BA.2 conv.')){
  print(srGroup)
  
  save_text <- strsplit(srGroup, "/")[[1]][1]
  if(srGroup == "AZ/BNT"){
    save_text <- "AZ_BNT"
  }
  newMap <- read.acmap(paste0("./som/bootstrapping/wo_",save_text,".ace"))
  srOutlineWidth(newMap) <- 1
  
  p <- procrustesMap(newMap, neut, sera = FALSE)
  
  title_text <- srGroup
  if(title_text == "mRNA1273/mRNA1273") {
    title_text <- "mRNA-1273/mRNA-1273"
  } else if (title_text == "WT conv.") {
    title_text <- "first wave conv."
  } else if (title_text == "BA.1 conv.") {
    title_text <- "BA.1 omicron conv."
  } else if (title_text == "BA.2 conv.") {
    title_text <- "BA.2 omicron conv."
  } else if (title_text == "alpha/alpha+E484K conv.") {
    title_text <- "alpha conv."
  }
  
  plot(p, xlim = xlim_no_zoom, ylim = ylim_no_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
       grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.3, plot_stress = TRUE)
  title(main = paste(title_text, 'sera', sep=' '), cex.main=0.7, line = 0.2)
  text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, labels[srGroup, ], cex = 1.2)
}

dev.off()
