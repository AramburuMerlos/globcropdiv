library(terra)
library(Recocrop)

# climatic data
tavg <- rast(Sys.glob("D:/WorldClim/2.1/wc5min/tavg/*.tif")) 
tmin <- rast(Sys.glob("D:/WorldClim/2.1/wc5min/tmin/*.tif")) 
prec <- rast(Sys.glob("D:/WorldClim/2.1/wc5min/prec/*.tif")) 

# ph
ph <- rast("D:/SoilGrids/phh2o/phh2o_0-15cm_mean_5min.tif")
ph <- ph/10

# crop mask
cm <- rast("D:/Dp_world/Monf_Cropland_Mask.tif")

# mask data
tavg <- mask(tavg, cm)
tmin <- mask(tmin, cm)
prec <- mask(prec, cm)
ph <- mask(ph, cm)

# crops to be modelled
dcrops <- read.csv(here::here("AuxData/Monfcrops.csv"), 
                  stringsAsFactors = F) 
crops <- unique(dcrops[, c("RID", "NAME", "SCIENTNAME")])
crops <- crops[order(crops$RID),]

dirn <- ("D:/Dp_world/CropSuit/")
dir.create(dirn)
dirn <- ("D:/Dp_world/CropSuit/Rainfed")
dir.create(dirn)
dirn <- ("D:/Dp_world/CropSuit/Rainfed/RID")
dir.create(dirn)


# check if all crop names or scientific names are found
isfound <- function(x) is.list(ecocropPars(x))
sink("NUL") # to avoid all messages
test <- all(((sapply(crops$NAME, isfound)|sapply(crops$SCIENTNAME, isfound))))
sink()
test

for(i in 1:nrow(crops)){
  if(is.list(ecocropPars(crops$NAME[i]))){
    crop <- ecocropPars(crops$NAME[i])
  } else {
    crop <- ecocropPars(crops$SCIENTNAME[i])
  }
  m <- ecocrop(crop)
  control(m, get_max=TRUE)
  fn <- paste0("ID", formatC(crops$RID[i], width = 4, flag = "0"))
  predict(m, tavg = tavg, ktmp = tmin, prec = prec, ph = ph, 
          filename = file.path(dirn, paste0(fn, ".tif")), 
          overwrite = T, 
          wopt = list(names = fn, filetype = "GTiff", 
                      gdal = c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))
}

