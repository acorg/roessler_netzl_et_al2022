library(Racmacs)

neut <- read.acmap('data/maps/map-OmicronI+II+III-thresholded-single_exposure-P1m1.ace')

xlim_no_zoom <- read.csv("./data/metadata/xlim_no_zoom.csv")$x
ylim_no_zoom <- read.csv("./data/metadata/ylim_no_zoom.csv")$x

# Resampling antigens and sera
neutBootAS <- bootstrapMap(
  neut,
  "resample",
  bootstrap_repeats = 1000,
  bootstrap_ags = TRUE,
  bootstrap_sr = TRUE,
  reoptimize = TRUE,
  optimizations_per_repeat = 100,
  ag_noise_sd = 0.7,
  titer_noise_sd = 0.7,
  options = list()
)

# Resampling antigens
neutBootA <- bootstrapMap(
  neut,
  "resample",
  bootstrap_repeats = 1000,
  bootstrap_ags = TRUE,
  bootstrap_sr = FALSE,
  reoptimize = TRUE,
  optimizations_per_repeat = 100,
  ag_noise_sd = 0.7,
  titer_noise_sd = 0.7,
  options = list()
)
# Resampling sera
neutBootS <- bootstrapMap(
  neut,
  "resample",
  bootstrap_repeats = 1000,
  bootstrap_ags = FALSE,
  bootstrap_sr = TRUE,
  reoptimize = TRUE,
  optimizations_per_repeat = 100,
  ag_noise_sd = 0.7,
  titer_noise_sd = 0.7,
  options = list()
)
# do the blobs
neutBootASBlobs <- bootstrapBlobs(neutBootAS, conf.level = 0.68, smoothing = 4, gridspacing = 0.05)
neutBootSBlobs <- bootstrapBlobs(neutBootS, conf.level = 0.68, smoothing = 4, gridspacing = 0.05)
neutBootABlobs <- bootstrapBlobs(neutBootA, conf.level = 0.68, smoothing = 4, gridspacing = 0.05)


# Plot figure
png("som/bootstrapping/resample-bootstrap.png", width = 9, height = 3, units = 'in', res=300, pointsize = 18)
layout(matrix(c(1, 2, 3), ncol=3))
par(oma=c(0, 0, 0, 0), mar=c(1, 0, 1, 0))
plot(neutBootASBlobs, xlim = xlim_no_zoom, ylim = ylim_no_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.7)
text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, 'A', cex = 1.4)
plot(neutBootABlobs, xlim = xlim_no_zoom, ylim = ylim_no_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.7)
text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, 'B', cex = 1.4)
plot(neutBootSBlobs, xlim = xlim_no_zoom, ylim = ylim_no_zoom, fill.alpha = 0.9, plot_labels = FALSE, outline.alpha = 0.9,
     grid.col = "#cfcfcf", grid.margin.col="#7d7d7d", cex=0.7)
text(xlim_no_zoom[1]+0.4, ylim_no_zoom[2]-0.4, 'C', cex = 1.4)

dev.off()


