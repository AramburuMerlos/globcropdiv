library(terra)

nwd <- nchar(getwd())
if(substr(getwd(),  nwd - 10, nwd) != "globcropdiv") stop("See 0000_wd.R")

dpath <- "D:/globcropdiv"
ppath <- "OutData"
Ddrive <- dir.exists(dpath) & dir.exists("D:/AQUASTAT")

# import rainfed suitability
rpath <- paste0(ifelse(Ddrive, dpath, ppath) ,"/EcocropSuit/byEcoID/Rainfed/RID")
rfn <- Sys.glob(file.path(rpath, "*.tif"))
rainfed <- rast(rfn)

# import irrigated suitability
ipath <- paste0(ifelse(Ddrive, dpath, ppath) ,"/EcocropSuit/byEcoID/Irrigated/RID")
ifn <- Sys.glob(file.path(ipath, "*.tif"))
irrigated <- rast(ifn)

# check if IDs match 
all.equal(sub(rpath, "", rfn), 
          sub(ipath, "", ifn))

ids <- sub(paste0(rpath,"/"),  "", sub(".tif", "", rfn))

# area equipped with irrigation as fraction of total cropland
iperc <- rast(paste0(ifelse(Ddrive, "D:/", "InData/"), "AQUASTAT/gmia_v5_aei_pct.asc"))
ifrac <- iperc/100

# function to compute weighted average
f <- function(rfed, irr, ifrac) rfed * (1-ifrac) + irr * ifrac

outpath <- file.path(ifelse(Ddrive, dpath, ppath), "EcocropSuit/byEcoID")

for(i in 1:length(ids)){
  lapp(c(rainfed[[i]], irrigated[[i]], ifrac), fun = f,
       filename = file.path(outpath, paste0(ids[[i]], ".tif")), 
       overwrite = T, 
       wopt = list(names = ids[[i]], filetype = "GTiff",
                   gdal = c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))
}
