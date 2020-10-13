library(terra)

nwd <- nchar(getwd())
if(substr(getwd(),  nwd - 10, nwd) != "globcropdiv") stop("See 0000_wd.R")

Ddrive <- dir.exists("D:/SoilGrids") 
if(Ddrive) {
  ph05 <- rast("D:/SoilGrids/phh2o/phh2o_0-5cm_mean.vrt")
  ph15 <- rast("D:/SoilGrids/phh2o/phh2o_5-15cm_mean.vrt")
} else {
  ph05 <- rast("InData/SoilGrids/phh2o/phh2o_0-5cm_mean.vrt")
  ph15 <- rast("InData/SoilGrids/phh2o/phh2o_5-15cm_mean.vrt")
}

ph015 <- c(ph05, ph15)
ph <- app(ph015, fun = mean, na.rm = T)

# model spatraster and file name
if(Ddrive){
  r <- rast("D:/globcropdiv/Monf_Cropland_Mask.tif") 
  outfn <- "D:/SoilGrids/phh2o/phh2o_0-15cm_mean_5min.tif"
} else {
  r <- rast("OutData/globcropdiv/Monf_Cropland_Mask.tif")
  outfn <- "InData/SoilGrids/phh2o/phh2o_0-15cm_mean_5min.tif"
}

# first aggregate to lower resolution to avoid NA 
# ph data at 250m resolution. Cropaland data at 5' arc (~5 and 8 km)
# factor will be 20 (5000/250)
ph_lr <- aggregate(ph, fact = 20, fun = "mean", na.rm = T)

# then project
project(ph_lr, r, overwrite = T, 
        filename = outfn, 
        wopt = list(progress = 1, memfrac = 0.8))
