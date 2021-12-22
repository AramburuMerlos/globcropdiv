library(magrittr)
library(data.table)
library(terra)
library(viridis)

if(system('hostname', TRUE) == "ESP-RH-9891"){
  setwd("D:/globcropdiv/")
} else if(system('hostname', TRUE) == "LAPTOP-ST129J47"){ 
  setwd("G:/My Drive/globcropdiv/")
} # else if { ... 


# import data ################

# country borders
countries <- geodata::world(resolution = 3, path = "InData/countries")


# tropics
w <- countries
# to use a faster buffering method (not lonlat)
crs(w) <- "+proj=utm +zone=1"
b <- buffer(w, 2)
trop <- vect(c("LINESTRING(-180 23, 180 23)", "LINESTRING(-180 -23, 180 -23)"))
trop <- erase(trop, b)
equ <- vect(c("LINESTRING(-180 0, 180 0)", "LINESTRING(-180 0, 180 0)"))
equ <- erase(equ, b)


# total cropland (ha) per cell
totcl <- rast("InData/TotalCropland.tif")
ar = ncol(totcl)/nrow(totcl)


# cropland by lat 
dlat <- fread("OutData/DfxLat.csv")
dlat[, Mha:= tot.area/1e6]


# plot lat -------------

tiff(filename = paste0('G:/My Drive/globcropdiv/Maps/Total_Cropland_lat.tif'), 
     units = "in",
     width = ncol(totcl)/600 * 1.2, 
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
mtext("a", line = -2, adj = 0.05, font = 2, cex = 1.5)


# lat
lat_tcl <- loess(Mha ~ lat, data = dlat, span = 0.1)

par(mar = c(3.5,2.5,2.8,1), xpd = NA)
plot(dlat$lat ~ dlat$Mha, 
     pch = 22, bg = "dark green", cex = 0.5, 
     xlab = "",
     ylab = "", 
     axes = F,
     xlim = c(0, 25),
     ylim = c(-82,95))
lines(dlat$lat ~ predict(lat_tcl),  col = "dark green", lwd = 3)
axis(side = 1, at = seq(0,25,5), cex.axis = 0.9)
axis(side = 2, at = seq(-90, 90, 30), cex.axis = 0.9)
mtext("Area (Mha)", side = 1, line = 2.2)
mtext("latidude", side = 2, line = 2.2)
mtext("b", line = -1.3, adj = 0.05, font = 2, cex = 1.5)
clip(-0.5, 25.05, -90, 90)
abline(h = c(-23, 23), lty = 3, lwd = 0.4)
abline(h = 0, lwd = 0.4)

dev.off()















# plot lat lon -------------

tiff(filename = paste0('G:/My Drive/globcropdiv/Maps/Total_Cropland.tif'), 
     units = "in",
     width = ncol(totcl)/600 * 1.1, 
     height = (ncol(totcl)/600)/ar * 1.5, 
     type = "cairo", res = 300, 
     compression = "zip")

layout(mat = matrix(c(1,3,2,0), ncol = 2), heights = c(3,2), widths = c(5,2))


# map
plot(totcl, type = "continuous", axes = FALSE, mar = c(2, 4, 2, 5),
     col = mako(256, direction = -1),
     maxcell = ncell(totcl))
plot(countries, lwd = 0.4, add = T)
mtext("a", line = -2, adj = 0.05, font = 2, cex = 1.5)


# lat
lat_tcl <- loess(Mha ~ lat, data = dlat, span = 0.1)

par(mar = c(3.5,2.5,2.5,1), xpd = NA)
plot(dlat$lat ~ dlat$Mha, 
     pch = 22, bg = "dark green", cex = 0.5, 
     xlab = "",
     ylab = "", 
     axes = F,
     xlim = c(0, 25),
     ylim = c(-82,95))
lines(dlat$lat ~ predict(lat_tcl),  col = "dark green", lwd = 3)
axis(side = 1, at = seq(0,25,5), cex.axis = 0.9)
axis(side = 2, at = seq(-90, 90, 30), cex.axis = 0.9)
mtext("Area (Mha)", side = 1, line = 2.2)
mtext("latidude", side = 2, line = 2.2)
mtext("b", line = -2, adj = 0.05, font = 2, cex = 1.5)

# lon 
lon_tcl <- loess(Mha ~ lon, data = dlon, span = 0.1)

par(mar = c(4,4,0,5), xpd = NA)
plot(dlon$lon, dlon$Mha, 
     pch = 22, bg = "dark green", cex = 0.5,
     xlab = "",
     ylab = "", 
     axes = F,
     xlim = c(-165, 165),
     ylim = c(0,17))
lines(dlon$lon, predict(lon_tcl),  col = "dark green", lwd = 3)
axis(side = 1, at = seq(-180,180,30), cex.axis = 0.9)
axis(side = 2, at = seq(0,20,5), cex.axis = 0.9)
mtext("longitude", side = 1, line = 2.2)
mtext("Area (Mha)", side = 2, line = 2.2)
mtext("c", line = -1, adj = 0.05, font = 2, cex = 1.5)

dev.off()







