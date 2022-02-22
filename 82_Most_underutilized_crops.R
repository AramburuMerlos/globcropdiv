library(terra)
library(RColorBrewer)
library(viridis)

#if(system('hostname', TRUE) == "ESP-RH-9891"){
setwd("D:/globcropdiv/")
#} else if(system('hostname', TRUE) == "LAPTOP-ST129J47"){ 
#  setwd("G:/My Drive/globcropdiv/")
#} # else if { ... 

# AVERAGE #########################################################

# upload -------

pcl <- rast("OutData/projected/CroplandProp.tif")

# mask low cropland cells
cl_mask <- pcl < 0.005

area <- rast("OutData/projected/CropPropGap.tif")
#area_eco <- rast("OutData/projected/CropPropGap_Eco.tif")
#area_sdm <- rast("OutData/projected/CropPropGap_SDM.tif")
crops <- names(area)

# identify -----

cell <- which(values(pcl)[,1] > 0.005)
m <- extract(area, cell, list = FALSE)

rs0 <- which(rowSums(m, na.rm = T) == 0) # remove these cells with no crop
cell <- cell[-rs0]
m <- m[-rs0,]
colnames(m) <- crops
# change to percentage
m <- m * 100

## 1st most underutilized  -----
# column with maximum proportion
col_1 <- max.col(m)
# most underutilized crop
crop_1 <- crops[col_1]
# proportion gap of most unutilized crop
prop_1 <- m[cbind(1:nrow(m), col_1)]


# reclassify minor crops ----------
cats <- read.csv("G:/My Drive/globcropdiv/AuxData/CropCategories2.csv")
#cats <- rbind(cats, data.frame(crop="none", crop.category="None", category=12))
cats <- cats[order(cats$category),]

cat_1 <- cats$category[match(crop_1, cats$crop)]


# create rasters -----
rc1 <- rp1 <- pcl
rv <- rep(NA, ncell(pcl))

rv[cell] <- cat_1
values(rc1) <- rv
writeRaster(
        rc1,
        overwrite = TRUE,
        filename = "OutData/underutil_cat.tif", 
        wopt = list(filetype = "GTiff", 
                    gdal=c("COMPRESS=Deflate","PREDICTOR=1", "ZLEVEL=6"))
)

rv[cell] <- prop_1
values(rp1) <- rv
writeRaster(
        rp1,
        overwrite = TRUE,
        filename = "OutData/underutil_prop.tif", 
        wopt = list(filetype = "GTiff", 
                    gdal=c("COMPRESS=Deflate","PREDICTOR=1", "ZLEVEL=6"))
)


# crop borders ----
box <- vect(c("LINESTRING(-150 72, 180 72)", "LINESTRING(-150 -57, 180 -57)"))
crs(box) <- "+proj=longlat +datum=WGS84 +no_defs"
box <- project(box, pcl)
ex <- ext(pcl)
ex[c(1,3,4)] <- ext(box)[c(1,3,4)]
ex[2] <- 15500000

rc1 <- crop(rc1, ex)
rp1 <- crop(rp1, ex)



# country borders ----
countries <- geodata::world(resolution = 3, path = "InData/countries")
countries <- project(countries, pcl)
countries <- crop(countries, ex)

# b <- buffer(countries, 1e5)
# trop <- vect(c("LINESTRING(-180 23, 180 23)", "LINESTRING(-180 -23, 180 -23)"))
# crs(trop) <- "+proj=longlat +datum=WGS84 +no_defs"
# trop <- project(trop, pcl)
# trop <- erase(trop, b)
# 
# equ <- vect(c("LINESTRING(-180 0, 180 0)", "LINESTRING(-180 0, 180 0)"))
# crs(equ) <- "+proj=longlat +datum=WGS84 +no_defs"
# equ <- project(equ, pcl)
# equ <- erase(equ, b)

# colors and brks -----

colors_cat <- brewer.pal(9, "Set1")

brks <- c(seq(0,25,5), 100)
colors_con <- mako(length(brks) + 1, begin = 0.3, direction = -1)
leg <- c(paste0(seq(0,20, 5), " - ", seq(5,25,5)), "> 25")


# plot x 1 ------
ar = ncol(rc1)/nrow(rc1)
tiff(filename = 'G:/My Drive/globcropdiv/Maps/under_crop.tif', 
     units = "in",
     width = ncol(rc1)/300 * 2, 
     height = (ncol(rc1)/300)/ar * 2, 
     type = "cairo", res = 300, 
     compression = "zip")
par(mfrow = c(2,2))


### cat 1 -----
plot(rc1, pax = list(lwd = 0, labels = FALSE), 
     type = "classes",
     col = colors_cat,
     maxcell = ncell(rc1), 
     mar = c(1,ar,1,ar),
     legend = FALSE
)
plot(countries, lwd = 0.4, add = T)
#plot(trop, col = "grey50", lty = 3, lwd = 0.5, add = T)
#plot(equ, col = "grey50", lwd = 0.4, add = T)

### legend cat -------
plot(box, axes = FALSE, col = "white", mar = c(0,0,0,0))
legend(legend = c(unique(cats$crop.category)),
       cex = 1.5, x = "left", fill = colors_cat,
       title = "Crop", title.adj = 0.5)

### prop 1 -----

plot(rp1, pax = list(lwd = 0, labels = FALSE), 
     type = "interval",
     col = colors_con,
     breaks = brks,
     maxcell = ncell(rc1), 
     mar = c(1,ar,1,ar),
     legend = FALSE
)
plot(countries, lwd = 0.4, add = T)
#plot(trop, col = "grey50", lty = 3, lwd = 0.5, add = T)
#plot(equ, col = "grey50", lwd = 0.4, add = T)

### legend prop ------

plot(box, axes = FALSE, col = "white", mar = c(0,0,0,0))
legend(legend = leg,
       cex = 1.5, x = "left", fill = colors_con,
       title = "Proportion gap (%)", title.adj = 0.5)

dev.off()

# SDM #########################################################

## upload -------

pcl <- rast("OutData/projected/CroplandProp.tif")

## mask low cropland cells
cl_mask <- pcl < 0.005

#area <- rast("OutData/projected/CropPropGap.tif")
#area_eco <- rast("OutData/projected/CropPropGap_Eco.tif")
area <- rast("OutData/projected/CropPropGap_SDM.tif")
crops <- names(area)

## identify -----

cell <- which(values(pcl)[,1] > 0.005)
m <- extract(area, cell, list = FALSE)

rs0 <- which(rowSums(m, na.rm = T) == 0) # remove these cells with no crop
cell <- cell[-rs0]
m <- m[-rs0,]
colnames(m) <- crops
## change to percentage
m <- m * 100

## 1st most underutilized  -----
# column with maximum proportion
col_1 <- max.col(m)
# most underutilized crop
crop_1 <- crops[col_1]
# proportion gap of most unutilized crop
prop_1 <- m[cbind(1:nrow(m), col_1)]


## reclassify minor crops ----------
cats <- read.csv("G:/My Drive/globcropdiv/AuxData/CropCategories2.csv")
#cats <- rbind(cats, data.frame(crop="none", crop.category="None", category=12))
cats <- cats[order(cats$category),]

cat_1 <- cats$category[match(crop_1, cats$crop)]


## create rasters -----
rc_sdm <- rp_sdm <- pcl
rv <- rep(NA, ncell(pcl))

rv[cell] <- cat_1
values(rc_sdm) <- rv
writeRaster(
        rc_sdm,
        overwrite = TRUE,
        filename = "OutData/underutil_cat_SDM.tif", 
        wopt = list(filetype = "GTiff", 
                    gdal=c("COMPRESS=Deflate","PREDICTOR=1", "ZLEVEL=6"))
)

rv[cell] <- prop_1
values(rp_sdm) <- rv
writeRaster(
        rp_sdm,
        overwrite = TRUE,
        filename = "OutData/underutil_prop_SDM.tif", 
        wopt = list(filetype = "GTiff", 
                    gdal=c("COMPRESS=Deflate","PREDICTOR=1", "ZLEVEL=6"))
)


## crop borders ----
box <- vect(c("LINESTRING(-150 72, 180 72)", "LINESTRING(-150 -57, 180 -57)"))
crs(box) <- "+proj=longlat +datum=WGS84 +no_defs"
box <- project(box, pcl)
ex <- ext(pcl)
ex[c(1,3,4)] <- ext(box)[c(1,3,4)]
ex[2] <- 15500000

rc_sdm <- crop(rc_sdm, ex)
rp_sdm <- crop(rp_sdm, ex)


# Eco #########################################################

## upload -------

pcl <- rast("OutData/projected/CroplandProp.tif")

## mask low cropland cells
cl_mask <- pcl < 0.005

#area <- rast("OutData/projected/CropPropGap.tif")
#area_eco <- rast("OutData/projected/CropPropGap_Eco.tif")
area <- rast("OutData/projected/CropPropGap_Eco.tif")
crops <- names(area)

## identify -----

cell <- which(values(pcl)[,1] > 0.005)
m <- extract(area, cell, list = FALSE)

rs0 <- which(rowSums(m, na.rm = T) == 0) # remove these cells with no crop
cell <- cell[-rs0]
m <- m[-rs0,]
colnames(m) <- crops
## change to percentage
m <- m * 100

## 1st most underutilized  -----
# column with maximum proportion
col_1 <- max.col(m)
# most underutilized crop
crop_1 <- crops[col_1]
# proportion gap of most unutilized crop
prop_1 <- m[cbind(1:nrow(m), col_1)]


## reclassify minor crops ----------
cat_1 <- cats$category[match(crop_1, cats$crop)]


## create rasters -----
rc_eco <- rp_eco <- pcl
rv <- rep(NA, ncell(pcl))

rv[cell] <- cat_1
values(rc_eco) <- rv
writeRaster(
        rc_eco,
        overwrite = TRUE,
        filename = "OutData/underutil_cat_Eco.tif", 
        wopt = list(filetype = "GTiff", 
                    gdal=c("COMPRESS=Deflate","PREDICTOR=1", "ZLEVEL=6"))
)

rv[cell] <- prop_1
values(rp_eco) <- rv
writeRaster(
        rp_eco,
        overwrite = TRUE,
        filename = "OutData/underutil_prop_Eco.tif", 
        wopt = list(filetype = "GTiff", 
                    gdal=c("COMPRESS=Deflate","PREDICTOR=1", "ZLEVEL=6"))
)


## crop borders ----
rc_eco <- crop(rc_eco, ex)
rp_eco <- crop(rp_eco, ex)




## country borders ----
countries <- geodata::world(resolution = 3, path = "InData/countries")
countries <- project(countries, pcl)
countries <- crop(countries, ex)




## plot x 2 ------
### colors and brks -----

colors_cat <- brewer.pal(9, "Set1")

brks <- c(seq(0,25,5), 100)
colors_con <- mako(length(brks), begin = 0.3, direction = -1)
leg <- c(paste0(seq(0,20, 5), " - ", seq(5,25,5)), "> 25")


ar = ncol(rc_eco)/nrow(rc_eco)
tiff(filename = 'G:/My Drive/globcropdiv/Maps/under_crop_Eco_SDM.tif', 
     unit = "in",
     width = ncol(rc1)/300 * 2, 
     height = (ncol(rc1)/300)/ar * 3, 
     type = "cairo", res = 300, 
     compression = "zip")
par(mfrow = c(3,2))


### cat eco -----
plot(rc_eco, pax = list(lwd = 0, labels = FALSE), 
     type = "classes",
     col = colors_cat,
     maxcell = ncell(rc_eco), 
     mar = c(1,ar,1,ar),
     legend = FALSE
)
plot(countries, lwd = 0.4, add = T)
#plot(trop, col = "grey50", lty = 3, lwd = 0.5, add = T)
#plot(equ, col = "grey50", lwd = 0.4, add = T)
mtext("a", line = -3, adj = .95, font = 2, cex = 3)

### prop eco -----

plot(rp_eco, pax = list(lwd = 0, labels = FALSE), 
     type = "interval",
     col = colors_con,
     breaks = brks,
     maxcell = ncell(rc_eco), 
     mar = c(1,ar,1,ar),
     legend = FALSE
)
plot(countries, lwd = 0.4, add = T)
#plot(trop, col = "grey50", lty = 3, lwd = 0.5, add = T)
#plot(equ, col = "grey50", lwd = 0.4, add = T)
mtext("b", line = -3, adj = .95, font = 2, cex = 3)


### cat sdm -----
plot(rc_sdm, pax = list(lwd = 0, labels = FALSE), 
     type = "classes",
     col = colors_cat,
     maxcell = ncell(rc_eco), 
     mar = c(1,ar,1,ar),
     legend = FALSE
)
plot(countries, lwd = 0.4, add = T)
#plot(trop, col = "grey50", lty = 3, lwd = 0.5, add = T)
#plot(equ, col = "grey50", lwd = 0.4, add = T)
mtext("c", line = -3, adj = .95, font = 2, cex = 3)


### prop sdm -----

plot(rp_sdm, pax = list(lwd = 0, labels = FALSE), 
     type = "interval",
     col = colors_con,
     breaks = brks,
     maxcell = ncell(rc_eco), 
     mar = c(1,ar,1,ar),
     legend = FALSE
)
plot(countries, lwd = 0.4, add = T)
#plot(trop, col = "grey50", lty = 3, lwd = 0.5, add = T)
#plot(equ, col = "grey50", lwd = 0.4, add = T)
mtext("d", line = -3, adj = .95, font = 2, cex = 3)


### legend cat -------
plot(box, axes = FALSE, col = "white", mar = c(0,0,0,0))
legend(legend = c(unique(cats$crop.category)),
       cex = 2.5, x = "top", fill = colors_cat, ncol = 2, 
       title = "Crop", title.adj = 0.5)

### legend prop ------

plot(box, axes = FALSE, col = "white", mar = c(0,0,0,0))
legend(legend = leg,
       cex = 2.5, x = "top", fill = colors_con, ncol = 2,
       title = "Proportion gap (%)", title.adj = 0.5)




dev.off()




















# OLD R ###############################################################

# plot x 2 ------
ar = ncol(rc1)/nrow(rc1)
tiff(filename = 'G:/My Drive/globcropdiv/Maps/unutil_crop_x2.tif', units = "in",
     width = ncol(rc1)/300 * 2, 
     height = (ncol(rc1)/300)/ar * 3, 
     type = "cairo", res = 300, 
     compression = "zip")
par(mfrow = c(3,2))


### cat 1 -----
plot(rc1, pax = list(lwd = 0, labels = FALSE), 
     type = "classes",
     col = colors_cat,
     maxcell = ncell(rc1), 
     mar = c(1,ar,1,ar),
     legend = FALSE
)
plot(countries, lwd = 0.4, add = T)
#plot(trop, col = "grey50", lty = 3, lwd = 0.5, add = T)
#plot(equ, col = "grey50", lwd = 0.4, add = T)


### prop 1 -----

plot(rp1, pax = list(lwd = 0, labels = FALSE), 
     type = "interval",
     col = colors_con,
     breaks = brks,
     maxcell = ncell(rc1), 
     mar = c(1,ar,1,ar),
     legend = FALSE
)
plot(countries, lwd = 0.4, add = T)
#plot(trop, col = "grey50", lty = 3, lwd = 0.5, add = T)
#plot(equ, col = "grey50", lwd = 0.4, add = T)


### cat 2 -----
plot(rc2, pax = list(lwd = 0, labels = FALSE), 
     type = "classes",
     col = colors_cat,
     maxcell = ncell(rc1), 
     mar = c(1,ar,1,ar),
     legend = FALSE
)
plot(countries, lwd = 0.4, add = T)
#plot(trop, col = "grey50", lty = 3, lwd = 0.5, add = T)
#plot(equ, col = "grey50", lwd = 0.4, add = T)


### prop 2 -----

plot(rp2, pax = list(lwd = 0, labels = FALSE), 
     type = "interval",
     col = colors_con,
     breaks = brks,
     maxcell = ncell(rc1), 
     mar = c(1,ar,1,ar),
     legend = FALSE
)
plot(countries, lwd = 0.4, add = T)
#plot(trop, col = "grey50", lty = 3, lwd = 0.5, add = T)
#plot(equ, col = "grey50", lwd = 0.4, add = T)


### legend -------
plot(box, axes = FALSE, col = "white")
legend(legend = c(unique(cats$crop.category)),
       cex = 3.3, x = "bottomleft", fill = colors_cat,
       title = "Crop", title.adj = 0.5, ncol = 2)

plot(box, axes = FALSE, col = "white")
legend(legend = leg,
       cex = 4, x = "bottomright", fill = colors_con,
       title = "Proportion gap (%)", title.adj = 0.5, ncol = 2)

dev.off()

