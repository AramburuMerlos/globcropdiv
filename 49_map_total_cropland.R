library(viridis)
library(terra)

r <- rast("InData/TotalCropland.tif")

ar = ncol(r)/nrow(r)

countries <- rnaturalearth::ne_download(scale = 10, type = "countries")

tiff(filename = paste0('Maps/Total_Cropland.tif'), 
     units = "in",
     width = ncol(r)/300, 
     height = (ncol(r)/300)/ar, 
     type = "cairo", res = 300, 
     compression = "zip")
par(mai = c(0,0,0,0), omi = c(0,0,0,0))
plot(r, type = "continuous", axes = FALSE,
     col = mako(256, direction = -1),
     maxcell = ncell(r), mar = c(2,0,2,6), 
     main = "Cropland (ha)")
plot(countries, lwd = 0.4, add = T)
dev.off()