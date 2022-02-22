library(terra)

if(system('hostname', TRUE) %in% c("ESP-RH-9891", "LAPTOP-ST129J47")){
  setwd("D:/globcropdiv/")
} # else if { ... 

world <- geodata::world(resolution = 1, path = "InData/countries")

r <- rast("InData/TotalCropland.tif")

rw <- rasterize(world, r, "GID_0")
writeRaster(rw,  
            filename = "InData/countries/rcountries.tif",
            overwrite = T, 
            wopt = list(names = "countries", filetype = "GTiff",
                        gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))

