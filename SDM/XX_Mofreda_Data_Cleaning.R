library(terra)

# to compute P/A estimations it might be safer to use only county level data. 
# here I create a version of the Monfreda dataset that keeps only county data
# Note that this reduces the total number of crops. 

nwd <- nchar(getwd())
if(substr(getwd(),  nwd - 10, nwd) != "globcropdiv") stop("See 0000_wd.R")

monpath <- if(dir.exists("D:/Monfreda")) "D:/Monfreda" else "InData/Monfreda"

# Crop Abundance Data and Data Quality ----------
abfn <- Sys.glob(file.path(monpath, "GeoTiff/*/*HarvestedAreaFraction.tif"))
abund <- rast(abfn)

dqfn <- Sys.glob(file.path(monpath, "GeoTiff/*/*DataQuality_HarvestedArea.tif"))
dq <- rast(dqfn)

crops <- sub(file.path(monpath, "GeoTiff/*"), "", 
             sub("/*_HarvestedAreaFraction.tif", "", abfn))
crops <- substr(crops, 1, ceiling(nchar(crops)/2)-1)

# Normalize all abundances to a max of 1
outpath <- file.path(monpath, "CountyLevel")
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


