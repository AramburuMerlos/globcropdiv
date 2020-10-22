library(terra)

# to compute P/A estimations it might be safer to use only county level data. 
# here I create a version of the Monfreda dataset that keeps only county data
# Note that this reduces the total number of crops. 

if(dir.exists("D:/Monfreda")){
  wcpath <- "D:/WorldClim/2.1/wc5min"
  sgpath <- "D:/SoilGrids"
  aqpath <- "D:/AQUASTAT"
  monpath <- "D:/Monfreda"
} else {
  wcpath <- "InData/WorldClim/2.1/wc5min"
  sgpath <- "InData/SoilGrids"
  aqpath <- "InData/AQUASTAT"
  monpath <- "InData/Monfreda"
}

# Crop Abundance Data and Data Quality ----------
abfn <- Sys.glob(file.path(monpath, "GeoTiff/*/*HarvestedAreaFraction.tif"))
abund <- rast(abfn)

dqfn <- Sys.glob(file.path(monpath, "GeoTiff/*/*DataQuality_HarvestedArea.tif"))
dq <- rast(dqfn)

crops <- sub(file.path(monpath, "GeoTiff/*"), "", 
             sub("/*_HarvestedAreaFraction.tif", "", abfn))
crops <- substr(crops, 1, ceiling(nchar(crops)/2)-1)

# Mask predictors NA values -----
preds <- c(rast(file.path(wcpath, "extra/AI.tif")),  
           rast(file.path(sgpath, "phh2o/phh2o_0-15cm_mean_5min.tif")), 
           rast(file.path(aqpath, "gmia_v5_aei_pct.asc"))  # irrigation
)

for(i in 1:nlyr(preds)){
  abund <- mask(abund, preds[[i]])
}


# Keep only County level data (dq == 1)
outpath <- file.path(monpath, "OnlyCounty")
dir.create(outpath)

for(i in 1:length(crops)){
  dqmax <- global(dq[[i]], "max", na.rm = T)[[1]]
  if(dqmax == 1){
    mask(abund[[i]], dq[[i]], inverse = T, maskvalue = 1,
         filename = file.path(outpath, paste0(crops[i], "_CountyLevel.tif")), 
         overwrite = T, 
         wopt = list(names = crops[i], filetype = "GTiff",
                     gdal = c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))
  }
}


