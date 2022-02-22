library(terra)
library(RColorBrewer)
library(viridis)

#if(system('hostname', TRUE) == "ESP-RH-9891"){
  setwd("D:/globcropdiv/")
#} else if(system('hostname', TRUE) == "LAPTOP-ST129J47"){ 
#  setwd("G:/My Drive/globcropdiv/")
#} # else if { ... 

# upload -------

pcl <- rast("OutData/projected/CroplandProp.tif")

# mask low cropland cells
cl_mask <- pcl < 0.005

area <- rast("OutData/projected/ActualCropProp.tif")
crops <- names(area)

# identify -----

cell <- which(values(pcl)[,1] > 0.005)
m <- extract(area, cell, list = FALSE)

rs0 <- which(rowSums(m, na.rm = T) == 0) # remove these cells with no crop
cell <- cell[-rs0]
m <- m[-rs0,]
colnames(m) <- crops

# need to replace NA for 0 to use max.col function
m <- replace(m, is.na(m), 0)

## 1st most dominant  -----
# column with maximum proportion
col_1 <- max.col(m)
# most dominant crop
crop_1 <- crops[col_1]
# proportion of most dominant crop
prop_1 <- m[cbind(1:nrow(m), col_1)]

# set most dominant crop proportion to 0 to get 2nd most dominant
m[cbind(1:nrow(m), col_1)] <- 0


## 2nd most dominant  -----
col_2 <- max.col(m)
crop_2 <- crops[col_2]
prop_2 <- m[cbind(1:nrow(m), col_2)]
# when there is no room for a dominant second crop
crop_2[prop_1 > 0.9] <- "none"

m[cbind(1:nrow(m), col_2)] <- 0

# reclassify minor crops ----------
cats <- read.csv("G:/My Drive/globcropdiv/AuxData/CropCategories2.csv")
cats <- cats[order(cats$category),]

cat_1 <- cats$category[match(crop_1, cats$crop)]
cat_2 <- cats$category[match(crop_2, cats$crop)]



# create rasters -----

rc1 <- rp1 <- rc2 <- rp2 <- pcl
rv <- rep(NA, ncell(pcl))

rv[cell] <- cat_1
values(rc1) <- rv
writeRaster(
  rc1,
  overwrite = TRUE,
  filename = "OutData/dominant_1_cat.tif", 
  wopt = list(filetype = "GTiff", 
              gdal=c("COMPRESS=Deflate","PREDICTOR=1", "ZLEVEL=6"))
)


rv[cell] <- prop_1
values(rp1) <- rv
writeRaster(
  rp1,
  overwrite = TRUE,
  filename = "OutData/dominant_1_prop.tif", 
  wopt = list(filetype = "GTiff", 
              gdal=c("COMPRESS=Deflate","PREDICTOR=1", "ZLEVEL=6"))
)

rv[cell] <- cat_2
values(rc2) <- rv
writeRaster(
  rc2, 
  filename = "OutData/dominant_2_cat.tif",
  overwrite = TRUE,
  wopt = list(filetype = "GTiff", 
              gdal=c("COMPRESS=Deflate","PREDICTOR=1", "ZLEVEL=6"))
)


rv[cell] <- prop_2
values(rp2) <- rv
writeRaster(
  rp2, 
  filename = "OutData/dominant_2_prop.tif",
  overwrite = TRUE,
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

rc2 <- crop(rc2, ex)
rp2 <- crop(rp2, ex)


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

# colors -----
# qualitative colors for crop categories
colors_cat <- brewer.pal(9, "Set1")
colors_con <- viridis(10)


# plot x 2 ------
ar = ncol(rc1)/nrow(rc1)
tiff(filename = 'G:/My Drive/globcropdiv/Maps/domin_crop_x2.tif', units = "in",
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
mtext("a", line = -3, adj = .95, font = 2, cex = 3)

### prop 1 -----

plot(rp1, pax = list(lwd = 0, labels = FALSE), 
     type = "interval",
     col = colors_con,
     breaks = seq(0,1,.1),
     maxcell = ncell(rc1), 
     mar = c(1,ar,1,ar),
     legend = FALSE
)
plot(countries, lwd = 0.4, add = T)
#plot(trop, col = "grey50", lty = 3, lwd = 0.5, add = T)
#plot(equ, col = "grey50", lwd = 0.4, add = T)
mtext("b", line = -3, adj = .95, font = 2, cex = 3)

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
mtext("c", line = -3, adj = .95, font = 2, cex = 3)

### prop 2 -----
rp12 <- rp1 + rp2 

plot(rp12, pax = list(lwd = 0, labels = FALSE), 
     type = "interval",
     col = colors_con,
     breaks = seq(0,1,.1), #seq(0,.5,.05),
     maxcell = ncell(rc1), 
     mar = c(1,ar,1,ar),
     legend = FALSE
)
plot(countries, lwd = 0.4, add = T)
#plot(trop, col = "grey50", lty = 3, lwd = 0.5, add = T)
#plot(equ, col = "grey50", lwd = 0.4, add = T)
mtext("d", line = -3, adj = .95, font = 2, cex = 3)


### legend -------
plot(box, axes = FALSE, col = "white")
legend(legend = c(unique(cats$crop.category)),
       cex = 3.3, x = "top", fill = colors_cat,
       title = "Crop", title.adj = 0.5, ncol = 2)

plot(box, axes = FALSE, col = "white")
# legend(legend = paste0(seq(0,90,10), " - ", seq(10,100,10)),
#        cex = 3, x = "topleft", fill = colors_con,
#        title = "Proportion 1 (%)", title.adj = 0.5, ncol = 2)
# 
legend(legend = paste0(seq(0,90,10), " - ", seq(10,100,10)),
       cex = 3.5, x = "top", fill = colors_con,
       title = "Proportion (%)", title.adj = 0.5, ncol = 2)

dev.off()



# singletons #########################################

ar = ncol(rc1)/nrow(rc1)
tiff(filename = 'G:/My Drive/globcropdiv/Maps/domin_crop_1.tif', units = "in",
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
     mar = c(1,ar,1,ar)/2,
     legend = FALSE
)
plot(countries, lwd = 0.4, add = T)
#plot(trop, col = "grey50", lty = 3, lwd = 0.5, add = T)
#plot(equ, col = "grey50", lwd = 0.4, add = T)
mtext("a", line = -2.5, adj = .95, font = 2, cex = 2)

plot(box, axes = FALSE, col = "white", mar = c(0,0,0,0))
legend(legend = c(unique(cats$crop.category)),
       cex = 1.5, x = "left", fill = colors_cat,
       title = "Crop", title.adj = 0.5, ncol = 1)


### prop 1 -----

plot(rp1, pax = list(lwd = 0, labels = FALSE), 
     type = "interval",
     col = colors_con,
     breaks = seq(0,1,.1),
     maxcell = ncell(rc1), 
     mar = c(1,ar,1,ar)/2,
     legend = FALSE
)
plot(countries, lwd = 0.4, add = T)
#plot(trop, col = "grey50", lty = 3, lwd = 0.5, add = T)
#plot(equ, col = "grey50", lwd = 0.4, add = T)
mtext("b", line = -2.5, adj = .95, font = 2, cex = 2)

### legend -------
plot(box, axes = FALSE, col = "white", mar = c(0,0,0,0))
# legend(legend = paste0(seq(0,90,10), " - ", seq(10,100,10)),
#        cex = 3, x = "topleft", fill = colors_con,
#        title = "Proportion 1 (%)", title.adj = 0.5, ncol = 2)
# 
legend(legend = paste0(seq(0,90,10), " - ", seq(10,100,10)),
       cex = 1.5, x = "left", fill = colors_con,
       title = "Proportion (%)", title.adj = 0.5, ncol = 1)


dev.off()



tiff(filename = 'G:/My Drive/globcropdiv/Maps/domin_crop_2.tif', units = "in",
     width = ncol(rc1)/300 * 2, 
     height = (ncol(rc1)/300)/ar * 2, 
     type = "cairo", res = 300, 
     compression = "zip")
par(mfrow = c(2,2))


### cat 2 -----
plot(rc2, pax = list(lwd = 0, labels = FALSE), 
     type = "classes",
     col = colors_cat,
     maxcell = ncell(rc1), 
     mar = c(1,ar,1,ar)/2,
     legend = FALSE
)
plot(countries, lwd = 0.4, add = T)
#plot(trop, col = "grey50", lty = 3, lwd = 0.5, add = T)
#plot(equ, col = "grey50", lwd = 0.4, add = T)
mtext("a", line = -2.5, adj = .95, font = 2, cex = 2)

### legend -------
plot(box, axes = FALSE, col = "white", mar = c(0,0,0,0))
legend(legend = c(unique(cats$crop.category)),
       cex = 1.5, x = "left", fill = colors_cat,
       title = "Crop", title.adj = 0.5)

### prop 2 -----
rp12 <- rp1 + rp2 

plot(rp12, pax = list(lwd = 0, labels = FALSE), 
     type = "interval",
     col = colors_con,
     breaks = seq(0,1,.1), #seq(0,.5,.05),
     maxcell = ncell(rc1), 
     mar = c(1,ar,1,ar)/2,
     legend = FALSE
)
plot(countries, lwd = 0.4, add = T)
#plot(trop, col = "grey50", lty = 3, lwd = 0.5, add = T)
#plot(equ, col = "grey50", lwd = 0.4, add = T)
mtext("b", line = -2.5, adj = .95, font = 2, cex = 2)


### legend -------
plot(box, axes = FALSE, col = "white", mar = c(0,0,0,0))
# legend(legend = paste0(seq(0,90,10), " - ", seq(10,100,10)),
#        cex = 3, x = "topleft", fill = colors_con,
#        title = "Proportion 1 (%)", title.adj = 0.5, ncol = 2)
# 
legend(legend = paste0(seq(0,90,10), " - ", seq(10,100,10)),
       cex = 1.5, x = "left", fill = colors_con,
       title = "Proportion (%)", title.adj = 0.5, ncol = 1)

dev.off()






# OLD R #####################################################



# plot x 3 ------
ar = ncol(rc1)/nrow(rc1)
tiff(filename = 'G:/My Drive/globcropdiv/Maps/domin_crop_x3.tif', units = "in",
     width = ncol(rc1)/300 * 2, 
     height = (ncol(rc1)/300)/ar * 4, 
     type = "cairo", res = 300, 
     compression = "zip")
par(mfrow = c(4,2))


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
     breaks = seq(0,1,.1),
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
     breaks = seq(0,1,.1),
     maxcell = ncell(rc1), 
     mar = c(1,ar,1,ar),
     legend = FALSE
)
plot(countries, lwd = 0.4, add = T)
#plot(trop, col = "grey50", lty = 3, lwd = 0.5, add = T)
#plot(equ, col = "grey50", lwd = 0.4, add = T)


### cat 3 -----
plot(rc3, pax = list(lwd = 0, labels = FALSE), 
     type = "classes",
     col = colors_cat,
     maxcell = ncell(rc1), 
     mar = c(1,ar,1,ar),
     legend = FALSE
)
plot(countries, lwd = 0.4, add = T)
#plot(trop, col = "grey50", lty = 3, lwd = 0.5, add = T)
#plot(equ, col = "grey50", lwd = 0.4, add = T)


### prop 3 -----

plot(rp3, pax = list(lwd = 0, labels = FALSE), 
     type = "interval",
     col = colors_con,
     breaks = seq(0,1,.1),
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
legend(legend = paste0(seq(0,90,10), " - ", seq(10,100,10)),
       cex = 4, x = "bottomright", fill = colors_con,
       title = "Proportion (%)", title.adj = 0.5, ncol = 2)

dev.off()

