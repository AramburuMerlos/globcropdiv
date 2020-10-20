library(terra)

nwd <- nchar(getwd())
if(substr(getwd(),  nwd - 10, nwd) != "globcropdiv") stop("See 0000_wd.R")

if(dir.exists("D:/globcropdiv")){
  wcpath <- "D:/WorldClim/2.1/wc5min"
  sgpath <- "D:/SoilGrids"
  gdpath <- "D:/globcropdiv"
  aqpath <- "D:/AQUASTAT"
  monpath <- "D:/Monfreda"
} else {
  wcpath <- "InData/WorldClim/2.1/wc5min"
  sgpath <- "InData/SoilGrids"
  aqpath <- "InData/AQUASTAT"
  gdpath <- "OutData"
  monpath <- "InData/Monfreda"
}

# Crop Abundance Data ----------
fn <- Sys.glob(file.path(monpath, "GeoTiff/*/*HarvestedAreaFraction.tif"))
crops <- sub(file.path(monpath, "GeoTiff/*"), "", 
             sub("/*_HarvestedAreaFraction.tif", "", fn))
crops <- substr(crops, 1, ceiling(nchar(crops)/2)-1)
abund <- rast(fn)

# cropland mask
cmask <- rast(file.path(gdpath, "Monf_Cropland_Mask.tif"))
# Predictors to mask NA values
bioc1 <- rast(Sys.glob(file.path(wcpath, "bioc/*.tif"))[1])
ph <- rast(file.path(sgpath, "phh2o/phh2o_0-15cm_mean_5min.tif"))
irr <- rast(file.path(aqpath, "gmia_v5_aei_pct.asc")) 

abund <- mask(abund, cmask)
abund <- mask(abund, bioc1)
abund <- mask(abund, ph)
abund <- mask(abund, irr)

# Normalize all abundances to a max of 1
outpath <- file.path(monpath, "Norm1")
dir.create(outpath)

for(i in 1:length(crops)){
  rmax <- global(abund[[i]], "max", na.rm = T)[[1]]
  if(rmax > 1){
    app(abund[[i]], fun = function(x) pmin(x, 1),
        filename = file.path(outpath, paste0(crops[i], "_norm1.tif")), 
        overwrite = T, 
        wopt = list(names = crops[i], filetype = "GTiff",
                    gdal = c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))
  } else if(rmax == 1){
    writeRaster(abund[[i]], 
        filename = file.path(outpath, paste0(crops[i], "_norm1.tif")), 
        overwrite = T, 
        wopt = list(names = crops[i], filetype = "GTiff",
                    gdal = c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))
  } else {
    app(abund[[i]], fun = function(x) x/rmax,
        filename = file.path(outpath, paste0(crops[i], "_norm1.tif")), 
        overwrite = T, 
        wopt = list(names = crops[i], filetype = "GTiff",
                    gdal = c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))
  }
}
  

