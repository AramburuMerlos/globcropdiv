library(magrittr)
library(terra)

if(!grepl("globcropdiv$", getwd())) warning("See 0000_wd.R")

# import rainfed suitability
rainfed <- "OutData/EcocropSuit/byEcoID/Rainfed/RID" %>% 
  file.path("*.tif") %>% Sys.glob %>% rast
  
# import irrigated suitability
irrigated <- "OutData/EcocropSuit/byEcoID/Irrigated/RID" %>% 
  file.path("*.tif") %>% Sys.glob %>% rast

# check if IDs match 
all.equal(names(irrigated), names(rainfed))

ids <- names(irrigated)

# area equipped with irrigation as fraction of total cropland
iperc <- rast("D:/AQUASTAT/gmia_v5_aei_pct.asc")
ifrac <- iperc/100

# function to compute weighted average
f <- function(rfed, irr, ifrac) rfed * (1-ifrac) + irr * ifrac

outpath <- file.path("OutData/EcocropSuit/byEcoID")

for(i in 1:length(ids)){
  lapp(c(rainfed[[i]], irrigated[[i]], ifrac), fun = f,
       filename = file.path(outpath, paste0(ids[[i]], ".tif")), 
       overwrite = T, 
       wopt = list(names = ids[[i]], filetype = "GTiff",
                   gdal = c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))
}
