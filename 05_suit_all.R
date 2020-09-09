library(terra)

# import rainfed suitability
rpath <- "D:/Dp_world/CropSuit/byEcoID/Rainfed/RID"
rfn <- Sys.glob(file.path(rpath, "*.tif"))
rainfed <- rast(rfn)

# import irrigated suitability
ipath <- "D:/Dp_world/CropSuit/byEcoID/Irrigated/RID/"
ifn <- Sys.glob(file.path(ipath, "*.tif"))
irrigated <- rast(ifn)

# check if IDs match 
all.equal(sub(rpath, "", rfn), 
          sub(ipath, "", ifn))

ids <- sub(paste0(rpath,"/"),  "", sub(".tif", "", rfn))

# area equipped with irrigation as fraction of total cropland
iperc <- rast("D:/AQUASTAT/gmia_v5_aei_pct.asc")
ifrac <- iperc/100

# function to compute weighted average
f <- function(rfed, irr, ifrac) rfed * (1-ifrac) + irr * ifrac

outpath <- "D:/Dp_world/CropSuit/byEcoID"

for(i in 1:length(ids)){
r <- lapp(c(rainfed[[i]], irrigated[[i]], ifrac), fun = f,
       filename = file.path(outpath, paste0(ids[[i]], ".tif")), 
       overwrite = T, 
       wopt = list(names = ids[[i]], filetype = "GTiff",
                   gdal = c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))
}
