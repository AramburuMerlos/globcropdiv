library(magrittr)
library(data.table)
library(terra)
library(RColorBrewer)

if(system('hostname', TRUE) == "ESP-RH-9891"){
  setwd("D:/globcropdiv/")
} else if(system('hostname', TRUE) == "LAPTOP-ST129J47"){ 
  setwd("G:/My Drive/globcropdiv/")
} # else if { ... 


# total cropland (ha) per cell
pcl <- rast("OutData/projected/CroplandProp.tif")

# mask low cropland cells
cl_mask <- pcl < 0.005

# crop borders
ex <- ext(pcl)
box <- vect(c("LINESTRING(-180 72, 180 72)", "LINESTRING(-180 -57, 180 -57)"))
crs(box) <- "+proj=longlat +datum=WGS84 +no_defs"
box <- project(box, pcl)
ex[1] <- -15500000
ex[3] <- ext(box)[3]
ex[4] <- ext(box)[4]


cl_mask <- crop(cl_mask, ex)

# country borders and tropic lines
countries <- geodata::world(resolution = 3, path = "InData/countries")
countries <- project(countries, pcl)


b <- buffer(countries, 1e5)
trop <- vect(c("LINESTRING(-180 23, 180 23)", "LINESTRING(-180 -23, 180 -23)"))
crs(trop) <- "+proj=longlat +datum=WGS84 +no_defs"
trop <- project(trop, pcl)
trop <- erase(trop, b)

equ <- vect(c("LINESTRING(-180 0, 180 0)", "LINESTRING(-180 0, 180 0)"))
crs(equ) <- "+proj=longlat +datum=WGS84 +no_defs"
equ <- project(equ, pcl)
equ <- erase(equ, b)


# colors -----
# qualitative colors for crop categories
colors_cat <- brewer.pal(9, "Set1")

# load -----
domin <- rast("OutData/dominant_1_cat.tif")
domin <- crop(domin, cl_mask)
domin <- mask(domin, cl_mask, maskvalues = 1)

under <- rast("OutData/underutil_cat_SDM.tif")
under <- crop(under, cl_mask)
under <- mask(under, cl_mask, maskvalues = 1)

cats <- read.csv("G:/My Drive/globcropdiv/AuxData/CropCategories2.csv")


# plot  ------
ar = ncol(domin)/nrow(domin)
tiff(filename = 'G:/My Drive/globcropdiv/Maps/domin_under.tif', units = "in",
     width = ncol(domin)/300, 
     height = (ncol(domin)/300)/ar * 2, 
     type = "cairo", res = 300, 
     compression = "zip")
par(mfrow = c(2,1))


### domin -----
plot(domin, pax = list(lwd = 0, labels = FALSE), 
     type = "classes",
     col = colors_cat,
     maxcell = ncell(domin), 
     mar = c(1,ar,1,ar)/2,
     legend = FALSE
)
plot(countries, lwd = 0.4, add = T)
plot(trop, col = "grey50", lty = 3, lwd = 0.5, add = T)
plot(equ, col = "grey50", lwd = 0.4, add = T)
mtext("a", line = -3, adj = .95, font = 2, cex = 3)


### under -----
plot(under, pax = list(lwd = 0, labels = FALSE), 
     type = "classes",
     col = colors_cat,
     maxcell = ncell(domin), 
     mar = c(1,ar,1,ar)/2,
     legend = FALSE
)
plot(countries, lwd = 0.4, add = T)
plot(trop, col = "grey50", lty = 3, lwd = 0.5, add = T)
plot(equ, col = "grey50", lwd = 0.4, add = T)
mtext("b", line = -3, adj = .95, font = 2, cex = 3)


### legend -------
legend(legend = c(unique(cats$crop.category))[c(1:3,5,4,6:9)],
       cex = 1.2, x = "bottomleft", fill = colors_cat[c(1:3,5,4,6:9)])#, 
       #title = "Crop", title.adj = 0.5, ncol = 1)

dev.off()


