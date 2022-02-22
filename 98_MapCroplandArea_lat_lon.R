library(magrittr)
library(data.table)
library(terra)
library(viridis)

if(system('hostname', TRUE) %in% c("ESP-RH-9891", "LAPTOP-ST129J47")){
  setwd("D:/globcropdiv/")
} # else if { ... 


# import data ################

# total cropland (ha) per cell
pcl <- rast("OutData/projected/CroplandProp.tif")
totcl <- pcl * prod(res(pcl))/1e4

# country borders
countries <- geodata::world(resolution = 3, path = "InData/countries")
countries <- project(countries, totcl)


b <- buffer(countries, 1e5)
trop <- vect(c("LINESTRING(-180 23, 180 23)", "LINESTRING(-180 -23, 180 -23)"))
crs(trop) <- "+proj=longlat +datum=WGS84 +no_defs"
trop <- project(trop, totcl)
trop <- erase(trop, b)

equ <- vect(c("LINESTRING(-180 0, 180 0)", "LINESTRING(-180 0, 180 0)"))
crs(equ) <- "+proj=longlat +datum=WGS84 +no_defs"
equ <- project(equ, totcl)
equ <- erase(equ, b)




ar = ncol(totcl)/nrow(totcl)

# cropland by lat 
dlat <- fread("OutData/DfxLat.csv")
dlat[, Mha:= tot.area/1e6]


# plot lat -------------

tiff(filename = paste0('G:/My Drive/globcropdiv/Maps/Total_Cropland_lat.tif'), 
     units = "in",
     width = ncol(totcl)/600 * 1.27, 
     height = (ncol(totcl)/600)/ar, 
     type = "cairo", res = 300, 
     compression = "zip")

layout(mat = matrix(c(1,2), ncol = 2), heights = 1, widths = c(5,2))


# map
plot(totcl, type = "continuous", mar = c(2, 4, 2, 5),
     pax = list(lwd = 0, labels = FALSE),
     col = mako(256, direction = -1),
     maxcell = ncell(totcl))
plot(countries, lwd = 0.4, add = T)
plot(trop, col = "grey50", lty = 3, lwd = 0.2, add = T)
plot(equ, col = "grey50", lwd = 0.2, add = T)
mtext("a", line = -1, adj = 0.05, font = 2, cex = 1.5)


# lat
lat_tcl <- loess(Mha ~ lat, data = dlat, span = 0.1)

par(mar = c(2.5,2.5,1.8,1), xpd = NA, mgp = c(3,0.5,0))
plot(dlat$lat ~ dlat$Mha, 
     pch = 22, bg = "dark green", cex = 0.5, 
     xlab = "",
     ylab = "", 
     axes = F,
     xlim = c(0, 25),
     ylim = c(-82,95))
lines(dlat$lat ~ predict(lat_tcl),  col = "dark green", lwd = 3)
axis(side = 1, at = seq(0,25,5), cex.axis = 0.8)
axis(side = 2, at = seq(-90, 90, 30), cex.axis = 0.8)
mtext("Area (Mha)", side = 1, line = 1.5)
mtext("latidude", side = 2, line = 1.5)
mtext("b", line = -2, adj = 0.05, font = 2, cex = 1.5)
clip(-0.5, 25.05, -90, 90)
abline(h = c(-23, 23), lty = 3, lwd = 0.4)
abline(h = 0, lwd = 0.4)

dev.off()

