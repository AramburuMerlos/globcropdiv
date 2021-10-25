library(magrittr)
library(data.table)
library(terra)
library(viridis)

if(system('hostname', TRUE) == "ESP-RH-9891"){
  setwd("D:/globcropdiv/")
} else if(system('hostname', TRUE) == "LAPTOP-ST129J47"){ 
  setwd("G:/My Drive/globcropdiv/")
} # else if { ... 

# country borders
countries <- geodata::world(resolution = 1, path = "InData/countries")
# total cropland (ha) per cell
totcl <- rast("InData/TotalCropland.tif")

# crop high latitudes
ex <- ext(countries)
ylim <- c(-55, 71)
ex[3] <- ylim[1]
ex[4] <- ylim[2]

b4 <- global(totcl, sum, na.rm = T)
totcl <- crop(totcl, ex)
af <- global(totcl, sum, na.rm = T)
all.equal(b4, af) # no cropland was left behind

# mask cells with less than 0.5% of cropland
area <- cellSize(totcl, unit = "ha")
pcl <- totcl/area
cl_mask <- pcl < 0.005

ar = ncol(totcl)/nrow(totcl)

# Total Cropland ----------------
totcl <- crop(totcl, ex)

tiff(filename = paste0('G:/My Drive/globcropdiv/Maps/Total_Cropland.tif'), 
     units = "in",
     width = ncol(totcl)/300, 
     height = (ncol(totcl)/300)/ar, 
     type = "cairo", res = 300, 
     compression = "zip")
par(mai = c(0,0,0,0), omi = c(0,0,0,0))
plot(totcl, type = "continuous", axes = FALSE,
     col = mako(256, direction = -1),
     maxcell = ncell(totcl))
plot(countries, lwd = 0.4, add = T)
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
plot(act_D, axes = FALSE, 
     type = "interval",
     col = mycolors,
     breaks = brks,
     #type = "continuous",
     #col = turbo(256, direction = -1),
     maxcell = ncell(act_D), 
     mar = c(0,0,0,0),
     plg = list(legend = c("1 - 1.5", "1.5 - 2.5", "2.5 - 5", "5 - 10", "> 10"),
                cex = 1.8, x = "bottomleft", 
                title = expression(italic(cD)),
                title.adj = 0.5))
plot(countries, lwd = 0.4, add = T)
dev.off()


## Potential #################################################


pot_D_eco <- rast("OutData/pot_D_eco.tif")
pot_D_eco <- crop(pot_D_eco, ex)
pot_D_eco <- mask(pot_D_eco, mask = cl_mask, maskvalues = 1)


pot_D_sdm <- rast("OutData/pot_D_sdm.tif")
pot_D_sdm <- crop(pot_D_sdm, ex)


tiff(filename = 'G:/My Drive/globcropdiv/Maps/pot_D.tif', units = "in",
     width = ncol(pot_D_sdm)/300, 
     height = (ncol(pot_D_sdm)/300)/ar * 2, 
     type = "cairo", res = 300, 
     compression = "zip")
par(mfrow = c(2,1))

### Ecocrop -------

brks <- seq(0, 175, length.out = 8)
mycols <- turbo(length(brks)-1, direction = -1, begin = 0.1, end = 0.9)

plot(pot_D_eco, axes = FALSE, 
     type = "interval",
     col = mycols,
     breaks = brks,
     maxcell = ncell(pot_D_eco), 
     mar = c(0,1,2,1),
     plg = list(cex = 1.8, x = "bottomleft", 
                title = expression("as-"*italic(pD)),
                title.adj = 0.3))
plot(countries, lwd = 0.4, add = T)
mtext("a", line = .1, adj = 0.95, cex = 2, font = 2)

### SDM -------

brks <- seq(50, 125, length.out = 7)
mycols <- turbo(length(brks)-1, direction = -1, begin = 0.1, end = 0.9)

plot(pot_D_sdm, axes = FALSE, 
     type = "interval",
     col = mycols,
     breaks = brks,
     maxcell = ncell(pot_D_sdm), 
     mar = c(0,1,2,1),
     plg = list(cex = 1.8,  x = "bottomleft", 
                title = expression("rs-"*italic(pD)),
                title.adj = 0.3))
plot(countries, lwd = 0.4, add = T)
mtext("b", line = .1, adj = 0.95, cex = 2, font = 2)

dev.off()








#### Species level ---------------
pot_D_eco_sp <- rast("OutData/pot_D_eco_sp.tif")
pot_D_eco_sp <- crop(pot_D_eco_sp, ex)
pot_D_eco_sp <- mask(pot_D_eco_sp, mask = cl_mask, maskvalues = 1)


brks <- c(brks, 200, 250, 350)
mycols <- turbo(length(brks)-1, direction = -1, begin = 0.1, end = 1)

tiff(filename = 'G:/My Drive/globcropdiv/Maps/pot_D_eco_sp.tif', units = "in",
     width = ncol(pot_D_eco_sp)/300, 
     height = (ncol(pot_D_eco_sp)/300)/ar, 
     type = "cairo", res = 300, 
     compression = "zip")
plot(pot_D_eco_sp, axes = FALSE, 
     type = "interval",
     col = turbo(length(brks)-1, direction = -1, begin = 0.1, end = 0.9),
     breaks = brks,
     maxcell = ncell(pot_D_eco_sp), 
     mar = c(0,0,0,0),
     plg = list(cex = 1.8, x = "left"))
plot(countries, lwd = 0.4, add = T)
dev.off()







## Attainable #################################################
brks <- c(1,5,10,15,20,25,30,50)

### SDM  -------------
att_D_sdm <- rast("OutData/att_D_sdm.tif")
att_D_sdm <- crop(att_D_sdm, ex)
att_D_sdm <- mask(att_D_sdm, cl_mask, maskvalues = 1)
mycols <- turbo(length(brks)-1, direction = -1, begin = 0.1, end = 1)


tiff(filename = 'G:/My Drive/globcropdiv/Maps/att_D.tif', units = "in",
     width = ncol(att_D_sdm)/300, 
     height = (ncol(att_D_sdm)/300)/ar * 2, 
     type = "cairo", res = 300, 
     compression = "zip")
par(mfrow = c(2,1))

plot(att_D_sdm, 
     axes = FALSE, 
     type = "interval", 
     col = mycols, 
     breaks = brks,
     maxcell = ncell(att_D_sdm), 
     mar = c(0,1,2,1),
     plg = list(legend = c("1 - 5", "5 - 10", "10 - 15", "15 - 20", 
                           "20 - 25", "25 - 30","> 30"),
                cex = 1.8, x = "bottomleft", 
                title = expression("rs-"*italic(aD)),
                title.adj = 0.3))
plot(countries, lwd = 0.4, add = T)
mtext("a", line = .2, cex = 2, font = 2, adj = 0.95)


### Ecocrop -------------
att_D_eco <- rast("OutData/att_D_eco.tif")
att_D_eco <- crop(att_D_eco, ex)
att_D_eco <- mask(att_D_eco, cl_mask, maskvalues = 1)


plot(att_D_eco, 
     axes = FALSE, 
     type = "interval", 
     col = mycols, 
     breaks = brks,
     maxcell = ncell(att_D_eco), 
     mar = c(0,1,2,1),
     plg = list(legend = c("1 - 5", "5 - 10", "10 - 15", "15 - 20", 
                           "20 - 25", "25 - 30","> 30"),
                cex = 1.8,  x = "bottomleft", 
                title = expression("as-"*italic(aD)),
                title.adj = 0.3))
plot(countries, lwd = 0.4, add = T)
mtext("b", line = .2, cex = 2, font = 2, adj = 0.95)
dev.off()










### Average -------------------------------
brks <- seq(0, 40, length.out = 6)
att_D_avg <- (att_D_eco + att_D_sdm)/2

att_D_avg <- crop(att_D_avg, ex)
att_D_avg <- mask(att_D_avg, cl_mask, maskvalues = 1)

tiff(filename = 'G:/My Drive/globcropdiv/Maps/att_D_avg.tif', units = "in",
     width = ncol(att_D_avg)/300, 
     height = (ncol(att_D_avg)/300)/ar, 
     type = "cairo", res = 300, 
     compression = "zip")
par(mai = c(0,0,0,0), omi = c(0,0,0,0))
plot(att_D_avg, 
     axes = FALSE, 
     type = "interval", 
     col = mycols, 
     breaks = brks,
     maxcell = ncell(att_D_avg), 
     mar = c(0,0,0,0),
     plg = list(legend = c("1 - 5", "5 - 10", "10 - 15", "15 - 20", 
                           "20 - 25", "25 - 30","> 30"),
                cex = 1.8,  x = "bottomleft", 
                title = expression(paste(phantom(0)^a,italic(D))),
                title.adj = 0.3))
plot(countries, lwd = 0.4, add = T)
dev.off()






## GAP ##################################################

brks <- c(-Inf, 30, 45, 60, 75, 90, 100)
mycols <- turbo(length(brks), begin = 0.1, end = 1)[-2]



### Ecocrop ---------
Dgap_eco <- (att_D_eco - act_D)/(att_D_eco) * 100

tiff(filename = 'G:/My Drive/globcropdiv/Maps/Dgap_eco.tif', units = "in",
     width = ncol(Dgap_eco)/300, 
     height = (ncol(Dgap_eco)/300)/ar, 
     type = "cairo", res = 300, 
     compression = "zip")
par(mai = c(0,0,0,0), omi = c(0,0,0,0))
plot(Dgap_eco, 
     axes = FALSE, 
     type = "interval", 
     col = mycols, 
     breaks = brks,
     maxcell = ncell(Dgap_eco), 
     mar = c(0,0,0,0),
     plg = list(cex = 1.8, x = "left"))
plot(countries, lwd = 0.45, add = T)
dev.off()

### SDM -------------
Dgap_sdm <-  (att_D_sdm - act_D)/(att_D_sdm) * 100


tiff(filename = 'G:/My Drive/globcropdiv/Maps/Dgap_sdm.tif', units = "in",
     width = ncol(Dgap_sdm)/300, 
     height = (ncol(Dgap_sdm)/300)/ar, 
     type = "cairo", res = 300, 
     compression = "zip")
par(mai = c(0,0,0,0), omi = c(0,0,0,0))
plot(Dgap_sdm, 
     axes = FALSE, 
     type = "interval", 
     col = mycols, 
     breaks = brks,
     maxcell = ncell(Dgap_sdm), 
     mar = c(0,0,0,0),
     plg = list(cex = 1.8, x = "left"))
plot(countries, lwd = 0.45, add = T)
dev.off()



### avg -------------
Dgap_avg <- (att_D_avg - act_D)/(att_D_avg) * 100


tiff(filename = 'G:/My Drive/globcropdiv/Maps/Dgap_avg.tif', units = "in",
     width = ncol(Dgap_avg)/300, 
     height = (ncol(Dgap_avg)/300)/ar, 
     type = "cairo", res = 300, 
     compression = "zip")
par(mai = c(0,0,0,0), omi = c(0,0,0,0))
plot(Dgap_avg, 
     axes = FALSE, 
     type = "interval", 
     col = mycols,
     breaks = brks,
     maxcell = ncell(Dgap_avg), 
     mar = c(0,0,0,0),
     plg = list(cex = 1.8, x = "bottomleft", 
                title = expression(italic(D)*g ~ "(%)"),
                title.adj = 0.3,
                legend = c("< 30", "30 - 45", "45 - 60", "60 - 75",
                           "75 - 90", "90 - 100")))
plot(countries, lwd = 0.45, add = T)
dev.off()

# DQI ###################################################
r <- totcl 

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
par(mai = c(0,0,0,0), omi = c(0,0,0,0))
plot(apple_dqi, type = "continuous", axes = FALSE,
     col = cividis(256, direction = -1),
     maxcell = ncell(r), mar = c(2,0,2,6), 
     plg = list(cex = 1.5))
plot(countries, lwd = 0.4, add = T)
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
plot(mango_dqi, type = "continuous", axes = FALSE,
     col = cividis(256, direction = -1),
     maxcell = ncell(r), mar = c(2,0,2,6), 
     plg = list(cex = 1.5))
plot(countries, lwd = 0.4, add = T)
mtext("b", line = .2, adj = .95, font = 2, cex = 2)
dev.off()




# Area and Suitability ###################################

## Actual  ----
dir_area <- 'G:/My Drive/globcropdiv/Maps/Area'
dir.create(dir_area, F, T)

dir_suit <- 'G:/My Drive/globcropdiv/Maps/Suit'
dir.create(dir_suit, F, T)

actual <- "InData/CropAbundance/*.tif" %>%  Sys.glob() %>%  rast()
crops <- names(actual)

d <- values(totcl) %>% 
  {data.table(cell = which(!is.na(.)), tcl = .[!is.na(.)])}
d[, (crops):= extract(actual, cell)]
d[, (crops):= lapply(.SD, `/`, tcl), .SDcols = crops]

# set 0 crop area as NA for plotting
for(j in crops) set(d, which(d[[j]] == 0), j, NA)

# empty raster
r <- rast(totcl)
v <- rep(NA_real_, ncell(r))

for(i in 1:length(crops)){
  v[d$cell] <- d[[crops[i]]]
  values(r) <- v
  
  tiff(filename = file.path(dir_area, paste0(crops[i], "_act.tif")), 
       units = "in",
       width = ncol(r)/300, 
       height = (ncol(r)/300)/ar, 
       type = "cairo", res = 300, 
       compression = "zip")
  par(mai = c(0,0,0,0), omi = c(0,0,0,0))
  plot(r, type = "continuous", axes = FALSE,
       col = mako(256, direction = -1),
       maxcell = ncell(r), mar = c(2,0,2,6), 
       main = paste(crops[i], "actual area (fraction)"))
  plot(countries, lwd = 0.4, add = T)
  dev.off()
}



## SDM Models ---------

### Data prep ---------
# Upload allocated SDM data
d <- fread("OutData/allocated_sdm.csv") 

acols <- names(d)[names(d) %like% '^a\\.']
scols <- gsub("^a", "s", acols)
crops <- gsub("^a\\.", "", acols)


# total allocated area
d[, tcl:= Reduce(`+`, .SD), .SDcols = acols]

# allocated area fraction
d[, (acols):= .SD/tcl, .SDcols = acols]

# set 0 crop area as NA for plotting
for(j in acols) set(d, which(d[[j]] == 0), j, NA)

# set 0 crop suit as NA for plotting
for(j in scols) set(d, which(d[[j]] == 0), j, NA)

# empty raster
r <- rast(totcl)
v <- rep(NA_real_, ncell(r))

### crop area ------------
for(i in 1:length(acols)){
  v[d$cell] <- d[[acols[i]]]
  values(r) <- v
  
  tiff(filename = file.path(dir_area, paste0(crops[i], "_sdm.tif")),  
       units = "in",
       width = ncol(r)/300, 
       height = (ncol(r)/300)/ar, 
       type = "cairo", res = 300, 
       compression = "zip")
  par(mai = c(0,0,0,0), omi = c(0,0,0,0))
  plot(r, type = "continuous", axes = FALSE,
       col = mako(256, direction = -1),
       maxcell = ncell(r), mar = c(2,0,2,6), 
       main = paste(crops[i], "allocated area (sdm)"))
  plot(countries, lwd = 0.4, add = T)
  dev.off()
}

### suitability -----------
for(i in 1:length(scols)){
  v[d$cell] <- d[[scols[i]]]
  values(r) <- v
  
  tiff(filename = file.path(dir_suit, paste0(crops[i], "_sdm.tif")),  
       units = "in",
       width = ncol(r)/300, 
       height = (ncol(r)/300)/ar, 
       type = "cairo", res = 300, 
       compression = "zip")
  par(mai = c(0,0,0,0), omi = c(0,0,0,0))
  plot(r, type = "continuous", axes = FALSE,
       col = mako(256, direction = -1),
       maxcell = ncell(r), mar = c(2,0,2,6), 
       main = paste(crops[i], "suitability (sdm)"))
  plot(countries, lwd = 0.4, add = T)
  dev.off()
}




## Ecocrop Model ----------

### Data prep ---------
# Upload allocated SDM data
d <- fread("OutData/allocated_eco.csv") 

acols <- names(d)[names(d) %like% '^a\\.']
scols <- gsub("^a", "s", acols)
crops <- gsub("^a\\.", "", acols)


# total allocated area
d[, tcl:= Reduce(`+`, .SD), .SDcols = acols]

# allocated area fraction
d[, (acols):= .SD/tcl, .SDcols = acols]

# set 0 crop area as NA for plotting
for(j in acols) set(d, which(d[[j]] == 0), j, NA)

# set 0 crop suit as NA for plotting
for(j in scols) set(d, which(d[[j]] == 0), j, NA)

# empty raster
r <- rast(totcl)
v <- rep(NA_real_, ncell(r))

### crop area ------------
for(i in 1:length(acols)){
  v[d$cell] <- d[[acols[i]]]
  values(r) <- v
  
  tiff(filename = file.path(dir_area, paste0(crops[i], "_eco.tif")),  
       units = "in",
       width = ncol(r)/300, 
       height = (ncol(r)/300)/ar, 
       type = "cairo", res = 300, 
       compression = "zip")
  par(mai = c(0,0,0,0), omi = c(0,0,0,0))
  plot(r, type = "continuous", axes = FALSE,
       col = mako(256, direction = -1),
       maxcell = ncell(r), mar = c(2,0,2,6), 
       main = paste(crops[i], "allocated area (eco)"))
  plot(countries, lwd = 0.4, add = T)
  dev.off()
}

### suitability -----------
for(i in 1:length(scols)){
  v[d$cell] <- d[[scols[i]]]
  values(r) <- v
  
  tiff(filename = file.path(dir_suit, paste0(crops[i], "_eco.tif")),  
       units = "in",
       width = ncol(r)/300, 
       height = (ncol(r)/300)/ar, 
       type = "cairo", res = 300, 
       compression = "zip")
  par(mai = c(0,0,0,0), omi = c(0,0,0,0))
  plot(r, type = "continuous", axes = FALSE,
       col = mako(256, direction = -1),
       maxcell = ncell(r), mar = c(2,0,2,6), 
       main = paste(crops[i], "suitability (eco)"))
  plot(countries, lwd = 0.4, add = T)
  dev.off()
}


