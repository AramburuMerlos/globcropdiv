library(magrittr)
library(data.table)
library(terra)
library(viridis)

if(system('hostname', TRUE) == "ESP-RH-9891"){
  setwd("D:/globcropdiv/")
} else if(system('hostname', TRUE) == "LAPTOP-ST129J47"){ 
  setwd("G:/My Drive/globcropdiv/")
} # else if { ... 


# total cropland (ha) per cell
pcl <- rast("OutData/projected/CroplandProp.tif")
totcl <- pcl * prod(res(pcl))/1e4

# mask low cropland cells
cl_mask <- pcl < 0.005
totcl <- mask(totcl, cl_mask, maskvalues = 1)

# crop borders
ex <- ext(pcl)
ylim <- vect(c("LINESTRING(-180 72, 180 72)", "LINESTRING(-180 -57, 180 -57)"))
crs(ylim) <- "+proj=longlat +datum=WGS84 +no_defs"
ylim <- project(ylim, totcl)
ex[1] <- -15500000
ex[3] <- ext(ylim)[3]
ex[4] <- ext(ylim)[4]


totcl <- crop(totcl, ex)
cl_mask <- crop(cl_mask, ex)

# country borders and tropic lines
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

# Total Cropland ----------------


tiff(filename = paste0('G:/My Drive/globcropdiv/Maps/Total_Cropland.tif'), 
     units = "in",
     width = ncol(totcl)/300, 
     height = (ncol(totcl)/300)/ar, 
     type = "cairo", res = 300, 
     compression = "zip")
plot(totcl, type = "continuous", axes = FALSE,
     col = mako(256, direction = -1),
     maxcell = ncell(totcl))
plot(countries, lwd = 0.4, add = T)
plot(trop, col = "grey50", lty = 3, lwd = 0.5, add = T)
plot(equ, col = "grey50", lwd = 0.4, add = T)
dev.off()


# Diversity ###########################################################

## Actual #############################################################
act_D <- rast("OutData/act_D.tif")
act_D <- crop(act_D, ex)
act_D <- mask(act_D, cl_mask, maskvalues = 1)

# pseudo - exponential breaks for better visualization
brks <- c(1, 1.5, 2.5, 5, 10, 50)
# select more colors and then select
mycolors <- turbo(10, direction = -1, begin = 0.1, end = 0.9)
mycolors <- mycolors[c(1,2,4,6,9)]

tiff(filename = 'G:/My Drive/globcropdiv/Maps/act_D.tif', units = "in",
     width = ncol(act_D)/300, 
     height = (ncol(act_D)/300)/ar, 
     type = "cairo", res = 300, 
     compression = "zip")
plot(act_D, pax = list(lwd = 0, labels = FALSE), 
     type = "interval",
     col = mycolors,
     breaks = brks,
     #type = "continuous",
     #col = turbo(256, direction = -1),
     maxcell = ncell(act_D), 
     mar = c(1,ar,1,ar),
     legend = FALSE
)
plot(countries, lwd = 0.4, add = T)
plot(trop, col = "grey50", lty = 3, lwd = 0.5, add = T)
plot(equ, col = "grey50", lwd = 0.4, add = T)
legend(legend = c("1 - 1.5", "1.5 - 2.5", "2.5 - 5", "5 - 10", "> 10"),
       cex = 1.3, x = "bottomleft", fill = mycolors,
       title = expression(italic(cD)),
       title.adj = 0.5)
dev.off()


## Potential #################################################
pot_D_eco <- rast("OutData/pot_D_eco.tif")
pot_D_eco <- crop(pot_D_eco, ex)
pot_D_eco <- mask(pot_D_eco, mask = cl_mask, maskvalues = 1)

pot_D_sdm <- rast("OutData/pot_D_sdm.tif")
pot_D_sdm <- crop(pot_D_sdm, ex)
pot_D_sdm <- mask(pot_D_sdm, mask = cl_mask, maskvalues = 1)


tiff(filename = 'G:/My Drive/globcropdiv/Maps/pot_D.tif', units = "in",
     width = ncol(pot_D_sdm)/300, 
     height = (ncol(pot_D_sdm)/300)/ar * 2, 
     type = "cairo", res = 300, 
     compression = "zip")
par(mfrow = c(2,1))


brks <- seq(0, 175, length.out = 8)
mycols <- turbo(length(brks)-1, direction = -1, begin = 0.1, end = 0.9)

plot(pot_D_eco, pax = list(lwd = 0, labels = FALSE), 
     type = "interval",
     col = mycols,
     breaks = brks,
     maxcell = ncell(pot_D_eco), 
     mar = c(1,ar,1,ar),
     legend = FALSE
)

plot(countries, lwd = 0.4, add = T)
plot(trop, col = "grey50", lty = 3, lwd = 0.5, add = T)
plot(equ, col = "grey50", lwd = 0.4, add = T)
legend(legend = paste(brks[-length(brks)], brks[-1], sep = " - "),
       fill = mycols,
       cex = 1.3, x = "bottomleft", 
       title = expression("as-"*italic(pD)),
       title.adj = 0.3)
mtext("a", line = -3, adj = 0.99, cex = 3, font = 2)

brks <- seq(50, 125, length.out = 7)
mycols <- turbo(length(brks)-1, direction = -1, begin = 0.1, end = 0.9)

plot(pot_D_sdm, pax = list(lwd = 0, labels = FALSE), 
     type = "interval",
     col = mycols,
     breaks = brks,
     maxcell = ncell(pot_D_sdm), 
     mar = c(1,ar,1,ar), 
     legend = FALSE
)

plot(countries, lwd = 0.4, add = T)
plot(trop, col = "grey50", lty = 3, lwd = 0.5, add = T)
plot(equ, col = "grey50", lwd = 0.4, add = T)
mtext("b", line = -3, adj = 0.99, cex = 3, font = 2)
legend(legend = paste(brks[-length(brks)], brks[-1], sep = " - "),
       fill = mycols,
       cex = 1.3,  x = "bottomleft", 
       title = expression("rs-"*italic(pD)),
       title.adj = 0.3)

dev.off()






## Attainable #################################################
brks <- c(1,5,10,15,20,25,30,50)

att_D_sdm <- rast("OutData/att_D_sdm.tif")
att_D_sdm <- crop(att_D_sdm, ex)
att_D_sdm <- mask(att_D_sdm, cl_mask, maskvalues = 1)

att_D_eco <- rast("OutData/att_D_eco.tif")
att_D_eco <- crop(att_D_eco, ex)
att_D_eco <- mask(att_D_eco, cl_mask, maskvalues = 1)


mycols <- turbo(length(brks)-1, direction = -1, begin = 0.1, end = 1)


tiff(filename = 'G:/My Drive/globcropdiv/Maps/att_D.tif', units = "in",
     width = ncol(att_D_sdm)/300, 
     height = (ncol(att_D_sdm)/300)/ar * 2, 
     type = "cairo", res = 300, 
     compression = "zip")
par(mfrow = c(2,1))

plot(att_D_sdm, 
     pax = list(lwd = 0, labels = FALSE), 
     type = "interval", 
     col = mycols, 
     breaks = brks,
     maxcell = ncell(att_D_sdm), 
     mar = c(1,ar,1,ar),
     legend = FALSE
)

plot(countries, lwd = 0.4, add = T)
plot(trop, col = "grey50", lty = 3, lwd = 0.5, add = T)
plot(equ, col = "grey50", lwd = 0.4, add = T)
mtext("a", line = -3, cex = 3, font = 2, adj = 0.99)
legend(legend = c("1 - 5", "5 - 10", "10 - 15", "15 - 20", 
                  "20 - 25", "25 - 30","> 30"),
       fill = mycols,
       cex = 1.3, x = "bottomleft",
       title = expression("rs-"*italic(aD)),
       title.adj = 0.3
)


plot(att_D_eco, 
     pax = list(lwd = 0, labels = FALSE), 
     type = "interval", 
     col = mycols, 
     breaks = brks,
     maxcell = ncell(att_D_eco), 
     mar = c(1,ar,1,ar), 
     legend = FALSE
)

plot(countries, lwd = 0.4, add = T)
plot(trop, col = "grey50", lty = 3, lwd = 0.5, add = T)
plot(equ, col = "grey50", lwd = 0.4, add = T)
mtext("b", line = -3, cex = 3, font = 2, adj = 0.99)
legend(legend = c("1 - 5", "5 - 10", "10 - 15", "15 - 20", 
                  "20 - 25", "25 - 30","> 30"),
       cex = 1.3,  x = "bottomleft", fill = mycols, 
       title = expression("as-"*italic(aD)),
       title.adj = 0.3
)
dev.off()







## GAP ##################################################

brks <- c(-Inf, 25, 50, 75, 100)
mycols <- turbo(length(brks)-1, begin = 0.1, end = .9)

### avg -------------
att_D_avg <- (att_D_eco + att_D_sdm)/2
Dgap_avg <- (att_D_avg - act_D)/(att_D_avg) * 100


tiff(filename = 'G:/My Drive/globcropdiv/Maps/Dgap_avg.tif', units = "in",
     width = ncol(Dgap_avg)/300, 
     height = (ncol(Dgap_avg)/300)/ar, 
     type = "cairo", res = 300, 
     compression = "zip")

plot(Dgap_avg, 
     pax = list(lwd = 0, labels = FALSE), 
     type = "interval", 
     col = mycols,
     breaks = brks,
     maxcell = ncell(Dgap_avg), 
     mar = c(1,ar,1,ar), 
     legend = FALSE
)
plot(countries, lwd = 0.45, add = T)
plot(trop, col = "grey50", lty = 3, lwd = 0.5, add = T)
plot(equ, col = "grey50", lwd = 0.4, add = T)
legend(cex = 1.3, x = "bottomleft", fill = mycols,
       title = expression(italic(D)*g ~ "(%)"),
       title.adj = 0.3,
       legend = c("< 25", "25 - 50", "50 - 75", "75 - 100")
) 
dev.off()





# DQI ###################################################
# unprojected total cropland
totcl <- rast("InData/TotalCropland.tif")
r <- totcl 

# unprojected countries
countries <- geodata::world(resolution = 3, path = "InData/countries")

apple_dqi <- rast("InData/CropAbundance/DataQualityIndex/apple_DQI.tif")
#apple <- rast("InData/CropAbundance/apple.tif")
#apple_dqi <- mask(apple_dqi, apple)
apple_dqi <- crop(apple_dqi, r)
apple_dqi <- 1 - apple_dqi/4

tiff(filename = file.path("G:/My Drive/globcropdiv/Maps/DQI_apple.tif"), 
     units = "in",
     width = ncol(r)/300, 
     height = (ncol(r)/300)/ar, 
     type = "cairo", res = 300, 
     compression = "zip")
par(mar = c(1,ar,1,ar), omi = c(0,0,0,0))
plot(apple_dqi, type = "continuous", pax = list(lwd = 0, labels = FALSE),
     col = cividis(256, direction = -1),
     maxcell = ncell(r), mar = c(2,0,2,6), 
     plg = list(cex = 1.5))
plot(countries, lwd = 0.4, add = T)
plot(trop, col = "grey50", lty = 3, lwd = 0.5, add = T)
plot(equ, col = "grey50", lwd = 0.4, add = T)
mtext("a", line = .2, adj = .95, font = 2, cex = 2)
dev.off()


mango_dqi <- rast("InData/CropAbundance/DataQualityIndex/mango_DQI.tif")
#mango <- rast("InData/CropAbundance/mango.tif")
#mango_dqi <- mask(mango_dqi, mango)
mango_dqi <- crop(mango_dqi, r)
mango_dqi <- 1 - mango_dqi/4

tiff(filename = file.path("G:/My Drive/globcropdiv/Maps/DQI_mango.tif"), 
     units = "in",
     width = ncol(r)/300, 
     height = (ncol(r)/300)/ar, 
     type = "cairo", res = 300, 
     compression = "zip")
par(mai = c(0,0,0,0), omi = c(0,0,0,0))
plot(mango_dqi, type = "continuous", pax = list(lwd = 0, labels = FALSE),
     col = cividis(256, direction = -1),
     maxcell = ncell(r), mar = c(2,0,2,6), 
     plg = list(cex = 1.5))
plot(countries, lwd = 0.4, add = T)
plot(trop, col = "grey50", lty = 3, lwd = 0.5, add = T)
plot(equ, col = "grey50", lwd = 0.4, add = T)
mtext("b", line = .2, adj = .95, font = 2, cex = 2)
dev.off()



# Diversity 2 ###############################################################


## Actual #############################################################
act_D2 <- rast("OutData/act_D2.tif")
act_D2 <- crop(act_D2, ex)
act_D2 <- mask(act_D2, cl_mask, maskvalues = 1)

# pseudo - exponential breaks for better visualization
brks <- c(1, 1.5, 2.5, 5, 10, 50)
# select more colors and then select
mycolors <- turbo(10, direction = -1, begin = 0.1, end = 0.9)
mycolors <- mycolors[c(1,2,4,6,9)]

tiff(filename = 'G:/My Drive/globcropdiv/Maps/act_D2.tif', units = "in",
     width = ncol(act_D2)/300, 
     height = (ncol(act_D2)/300)/ar, 
     type = "cairo", res = 300, 
     compression = "zip")
plot(act_D2, pax = list(lwd = 0, labels = FALSE), 
     type = "interval",
     col = mycolors,
     breaks = brks,
     #type = "continuous",
     #col = turbo(256, direction = -1),
     maxcell = ncell(act_D2), 
     mar = c(1,ar,1,ar),
     legend = FALSE
)
plot(countries, lwd = 0.4, add = T)
plot(trop, col = "grey50", lty = 3, lwd = 0.5, add = T)
plot(equ, col = "grey50", lwd = 0.4, add = T)
legend(legend = c("1 - 1.5", "1.5 - 2.5", "2.5 - 5", "5 - 10", "> 10"),
       cex = 1.3, x = "bottomleft", fill = mycolors,
       title = expression(italic(current)~phantom(0)^2*italic(D)),
       title.adj = 0.5)
dev.off()


## Potential #################################################
pot_D2_eco <- rast("OutData/pot_D2_eco.tif")
pot_D2_eco <- crop(pot_D2_eco, ex)
pot_D2_eco <- mask(pot_D2_eco, mask = cl_mask, maskvalues = 1)

pot_D2_sdm <- rast("OutData/pot_D2_sdm.tif")
pot_D2_sdm <- crop(pot_D2_sdm, ex)
pot_D2_sdm <- mask(pot_D2_sdm, mask = cl_mask, maskvalues = 1)


tiff(filename = 'G:/My Drive/globcropdiv/Maps/pot_D2.tif', units = "in",
     width = ncol(pot_D2_sdm)/300, 
     height = (ncol(pot_D2_sdm)/300)/ar * 2, 
     type = "cairo", res = 300, 
     compression = "zip")
par(mfrow = c(2,1))


brks <- seq(0, 175, length.out = 8)
mycols <- turbo(length(brks)-1, direction = -1, begin = 0.1, end = 0.9)

plot(pot_D2_eco, pax = list(lwd = 0, labels = FALSE), 
     type = "interval",
     col = mycols,
     breaks = brks,
     maxcell = ncell(pot_D2_eco), 
     mar = c(1,ar,1,ar),
     legend = FALSE
)

plot(countries, lwd = 0.4, add = T)
plot(trop, col = "grey50", lty = 3, lwd = 0.5, add = T)
plot(equ, col = "grey50", lwd = 0.4, add = T)
legend(legend = paste(brks[-length(brks)], brks[-1], sep = " - "),
       fill = mycols,
       cex = 1.3, x = "bottomleft", 
       title = expression("as-"*italic(p)*phantom(0)^2*italic(D)),
       title.adj = 0.3)
mtext("a", line = -3, adj = 0.99, cex = 3, font = 2)

brks <- seq(40, 120, length.out = 9)
mycols <- turbo(length(brks)-1, direction = -1, begin = 0.1, end = 0.9)

plot(pot_D2_sdm, pax = list(lwd = 0, labels = FALSE), 
     type = "interval",
     col = mycols,
     breaks = brks,
     maxcell = ncell(pot_D2_sdm), 
     mar = c(1,ar,1,ar), 
     legend = FALSE
)

plot(countries, lwd = 0.4, add = T)
plot(trop, col = "grey50", lty = 3, lwd = 0.5, add = T)
plot(equ, col = "grey50", lwd = 0.4, add = T)
mtext("b", line = -3, adj = 0.99, cex = 3, font = 2)
legend(legend = paste(brks[-length(brks)], brks[-1], sep = " - "),
       fill = mycols,
       cex = 1.3,  x = "bottomleft", 
       title = expression("rs-"*italic(p)*phantom(0)^2*italic(D)),
       title.adj = 0.3)

dev.off()






## Attainable #################################################
brks <- c(1,4,8,12,16,20,30)

att_D2_sdm <- rast("OutData/att_D2_sdm.tif")
att_D2_sdm <- crop(att_D2_sdm, ex)
att_D2_sdm <- mask(att_D2_sdm, cl_mask, maskvalues = 1)

att_D2_eco <- rast("OutData/att_D2_eco.tif")
att_D2_eco <- crop(att_D2_eco, ex)
att_D2_eco <- mask(att_D2_eco, cl_mask, maskvalues = 1)


mycols <- turbo(length(brks)-1, direction = -1, begin = 0.1, end = 1)


tiff(filename = 'G:/My Drive/globcropdiv/Maps/att_D2.tif', units = "in",
     width = ncol(att_D2_sdm)/300, 
     height = (ncol(att_D2_sdm)/300)/ar * 2, 
     type = "cairo", res = 300, 
     compression = "zip")
par(mfrow = c(2,1))

plot(att_D2_sdm, 
     pax = list(lwd = 0, labels = FALSE), 
     type = "interval", 
     col = mycols, 
     breaks = brks,
     maxcell = ncell(att_D2_sdm), 
     mar = c(1,ar,1,ar),
     legend = FALSE
)

plot(countries, lwd = 0.4, add = T)
plot(trop, col = "grey50", lty = 3, lwd = 0.5, add = T)
plot(equ, col = "grey50", lwd = 0.4, add = T)
mtext("a", line = -3, cex = 3, font = 2, adj = 0.99)
legend(legend = c("1 - 4", "4 - 8", "8 - 12", "12 - 16", 
                  "16 - 20", "> 20"),
       fill = mycols,
       cex = 1.3, x = "bottomleft",
       title = expression("rs-"*italic(a)*phantom(0)^2*italic(D)),
       title.adj = 0.3
)


plot(att_D2_eco, 
     pax = list(lwd = 0, labels = FALSE), 
     type = "interval", 
     col = mycols, 
     breaks = brks,
     maxcell = ncell(att_D2_eco), 
     mar = c(1,ar,1,ar), 
     legend = FALSE
)

plot(countries, lwd = 0.4, add = T)
plot(trop, col = "grey50", lty = 3, lwd = 0.5, add = T)
plot(equ, col = "grey50", lwd = 0.4, add = T)
mtext("b", line = -3, cex = 3, font = 2, adj = 0.99)
legend(legend = c("1 - 4", "4 - 8", "8 - 12", "12 - 16", 
                  "16 - 20", "> 20"),
       cex = 1.3,  x = "bottomleft", fill = mycols, 
       title = expression("as-"*italic(a)*phantom(0)^2*italic(D)),
       title.adj = 0.3
)
dev.off()







## GAP ##################################################

brks <- c(-Inf, 25, 50, 75, 100)
mycols <- turbo(length(brks)-1, begin = 0.1, end = .9)

### avg -------------
att_D2_avg <- (att_D2_eco + att_D2_sdm)/2
D2gap_avg <- (att_D2_avg - act_D2)/(att_D2_avg) * 100


tiff(filename = 'G:/My Drive/globcropdiv/Maps/D2gap_avg.tif', units = "in",
     width = ncol(D2gap_avg)/300, 
     height = (ncol(D2gap_avg)/300)/ar, 
     type = "cairo", res = 300, 
     compression = "zip")

plot(D2gap_avg, 
     pax = list(lwd = 0, labels = FALSE), 
     type = "interval", 
     col = mycols,
     breaks = brks,
     maxcell = ncell(D2gap_avg), 
     mar = c(1,ar,1,ar), 
     legend = FALSE
)
plot(countries, lwd = 0.45, add = T)
plot(trop, col = "grey50", lty = 3, lwd = 0.5, add = T)
plot(equ, col = "grey50", lwd = 0.4, add = T)
legend(cex = 1.3, x = "bottomleft", fill = mycols,
       title = expression(phantom(0)^2*italic(D)*g ~ "(%)"),
       title.adj = 0.3,
       legend = c("< 25", "25 - 50", "50 - 75", "75 - 100")
) 
dev.off()
