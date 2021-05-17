library(terra)
library(data.table)

if(system('hostname', TRUE) == "ESP-RH-9891"){
  setwd("D:/globcropdiv/")
} else if(system('hostname', TRUE) == "LAPTOP-ST129J47"){ 
  setwd("G:/My Drive/globcropdiv/")
} # else if { ... 


# pH ----

spath <- "InData/SoilGrids/phh2o/"

## High Resolution Average -----
ph05 <- rast(file.path(spath, "phh2o_0-5cm_mean.vrt"))
ph15 <- rast(file.path(spath, "phh2o_5-15cm_mean.vrt"))

ph015 <- c(ph05, ph15)
app(ph015, fun = mean, na.rm = T, 
    filename = file.path(spath, "pH0-15HR.tif"), 
    overwrite = T, 
    wopt = list(names = "pH0-15HR", filetype = "GTiff",
                gdal = c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))


## Lower the resolution -------
ph <- rast(file.path(spath, "pH0-15HR.tif"))

# first aggregate to lower resolution to avoid getting many NA 
# ph data is at 250m resolution. Cropland data is at 5' arc (~5 and 8 km)
# factor will be 20 (5000/250)

ph_lr <- aggregate(ph, fact = 20, fun = "mean", na.rm = T)

# then project

# raster to use as model for projection
r <- rast(Sys.glob("InData/WorldClim/2.1/wc5min/tavg/*.tif")[1])

outfn <- file.path(spath, "phh2o_0-15cm_mean_5min.tif")

project(ph_lr, r, overwrite = T, filename = outfn, 
        wopt = list(names = "pH0-15", filetype = "GTiff",
                    gdal = c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))
