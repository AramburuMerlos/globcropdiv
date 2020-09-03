library(terra)
ph05 <- rast("D:/SoilGrids/phh2o/phh2o_0-5cm_mean.vrt")
ph15 <- rast("D:/SoilGrids/phh2o/phh2o_5-15cm_mean.vrt")
ph015 <- c(ph05, ph15)
ph <- app(ph015, fun = mean)

# model spatraster 
r <- rast("D:/Dp_world/Monf_Cropland_Mask.tif") 

project(ph, r, overwrite = T, 
        filename = "D:/SoilGrids/phh2o/phh2o_0-15cm_mean_5min.tif", 
        wopt = list(progress = 1, memfrac = 0.8))
