library(terra)
library(Recocrop)

# climatic data
tavg <- rast(Sys.glob("D:/WorldClim/2.1/wc5min/tavg/*.tif")) 
tmin <- rast(Sys.glob("D:/WorldClim/2.1/wc5min/tmin/*.tif")) 

# ph
ph <- rast("D:/SoilGrids/phh2o/phh2o_0-15cm_mean_5min.tif")
ph <- ph/10
ph <- resample(ph, tavg[[1]]) # idk why this is necessary, but it is.

# irrigated crops mask
# AQUASTAT
ic <- rast("D:/AQUASTAT/gmia_v5_aei_pct.asc")
# consider only land with more than 20% of area irrigated
icm <- app(ic, function(x) x >= 20, nodes = 6)

# mask data
tavg <- mask(tavg, icm, maskvalue = 0)
tmin <- mask(tmin, icm, maskvalue = 0)
ph <- mask(ph, icm, maskvalue = 0)

# crops to be modelled
dcrops <- read.csv(here::here("AuxData/Monfcrops.csv"), 
                  stringsAsFactors = F) 
crops <- unique(dcrops[, c("RID", "NAME", "SCIENTNAME")])
crops <- crops[order(crops$RID),]

dirn <- ("D:/Dp_world/CropSuit/Irrigated")
dir.create(dirn)
dirn <- ("D:/Dp_world/CropSuit/Irrigated/RID")
dir.create(dirn)


# check if all crop names or scientific names are found
isfound <- function(x) is.list(ecocropPars(x))
sink("NUL") # to avoid all messages
test <- all(((sapply(crops$NAME, isfound)|sapply(crops$SCIENTNAME, isfound))))
sink()
test

sink("NUL")
for(i in 1:nrow(crops)){
  if(is.list(ecocropPars(crops$NAME[i]))){
    crop <- ecocropPars(crops$NAME[i])
  } else {
    crop <- ecocropPars(crops$SCIENTNAME[i])
  }
  m <- ecocrop(crop)
  control(m, get_max=TRUE)
  fn <- paste0("ID", formatC(crops$RID[i], width = 4, flag = "0"))
  predict(m, tavg = tavg, ktmp = tmin, ph = ph, 
          filename = file.path(dirn, paste0(fn, ".tif")), 
          overwrite = T, 
          wopt = list(names = fn, filetype = "GTiff",
          gdal = "COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))
}
sink()
