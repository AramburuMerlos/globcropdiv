library(terra)

nwd <- nchar(getwd())
if(substr(getwd(),  nwd - 10, nwd) != "globcropdiv") stop("See 0000_wd.R")

Ddrive <- dir.exists("D:/SoilGrids")
spath <- if(Ddrive) "D:/SoilGrids/phh2o/" else "InData/SoilGrids/phh2o/"
# -----
ph05 <- rast(file.path(spath, "phh2o_0-5cm_mean.vrt"))
ph15 <- rast(file.path(spath, "phh2o_5-15cm_mean.vrt"))

ph015 <- c(ph05, ph15)
app(ph015, fun = mean, na.rm = T, 
    filename = file.path(spath, "pH0-15HR.tif"), 
    overwrite = T, 
    wopt = list(names = "pH0-15HR", filetype = "GTiff",
                gdal = c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))
# -------
ph <- rast(file.path(spath, "pH0-15HR.tif"))

# model spatraster and file name
r <- rast(paste0(ifelse(Ddrive, "D:/", "OutData/"), 
                 "globcropdiv/Monf_Cropland_Mask.tif")) 
outfn <- file.path(spath, "phh2o_0-15cm_mean_5min.tif")

# first aggregate to lower resolution to avoid NA 
# ph data at 250m resolution. Cropaland data at 5' arc (~5 and 8 km)
# factor will be 20 (5000/250)
ph_lr <- aggregate(ph, fact = 20, fun = "mean", na.rm = T)
names(ph_lr) <- "pH0_15"
# then project
project(ph_lr, r, overwrite = T, 
        filename = outfn, 
        wopt = list(progress = 1, memfrac = 0.8, names = "pH0_15"))
