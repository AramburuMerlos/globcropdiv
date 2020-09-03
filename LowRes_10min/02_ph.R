library(terra)
ph05 <- rast("D:/SoilGrids/phh2o/phh2o_0-5cm_mean.vrt")
# model spatraster (climatic data)
r <- rast(Sys.glob("D:/WorldClim/2.1/wc10min/tavg/*.tif")[1]) 
project(ph05, r, overwrite = T, 
        filename = "D:/SoilGrids/phh2o/phh20_0-5cm_10min.tif", 
        wopt = list(progress = 1, memfrac = 0.6))
