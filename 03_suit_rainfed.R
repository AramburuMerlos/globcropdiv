library(terra)
library(Recocrop)

# climatic data
tavg <- rast(Sys.glob("D:/WorldClim/2.1/wc5min/tavg/*.tif")) 
tmin <- rast(Sys.glob("D:/WorldClim/2.1/wc5min/tmin/*.tif")) 
prec <- rast(Sys.glob("D:/WorldClim/2.1/wc5min/prec/*.tif")) 
annprec <- app(prec, fun = sum, na.rm = T)

# ph
ph <- rast("D:/SoilGrids/phh2o/phh2o_0-15cm_mean_5min.tif")
ph <- ph/10

# crop mask
cm <- rast("D:/Dp_world/Monf_Cropland_Mask.tif")

# mask data
tavg <- mask(tavg, cm)
tmin <- mask(tmin, cm)
prec <- mask(prec, cm)
annprec <- mask(annprec, cm)
ph <- mask(ph, cm)

# crops to be modelled
dcrops <- read.csv(here::here("AuxData/Monfcrops.csv"), 
                  stringsAsFactors = F) 
crops <- unique(dcrops[, c("RID", "NAME", "SCIENTNAME")])
crops <- crops[order(crops$RID),]

dirn <- ("D:/Dp_world/CropSuit/")
dir.create(dirn)
dirn <- ("D:/Dp_world/CropSuit/byEcoID")
dir.create(dirn)
dirn <- ("D:/Dp_world/CropSuit/byEcoID/Rainfed")
dir.create(dirn)
dirn <- ("D:/Dp_world/CropSuit/byEcoID/Rainfed/RID")
dir.create(dirn)

# ecocrop data base
fname <- system.file("parameters/ecocrop.rds", package="Recocrop")
ecodb <- readRDS(fname)
ecodb$NAME <- gsub("\\*", "", ecodb$NAME)
ecodb$NAME <- trimws(ecodb$NAME)
ecodb$nms <- paste(ecodb[,1], "; ", ecodb[,3], paste0(" (", ecodb[,2], ")"), sep="") 
raincols <- match(c("RMIN","ROPMN","ROPMX","RMAX"), names(ecodb))

# check if all Monf crops are found (by name or scientific name)
# 1. in ecocropPars function
isfound <- function(x) is.list(ecocropPars(x))
sink("NUL") # to avoid all messages
test <- all(((sapply(crops$NAME, isfound)|sapply(crops$SCIENTNAME, isfound))))
sink()
test

# 2. in EcoCrop Data Base
test2 <- vector(length = nrow(crops))
sink("NUL")
for (i in 1:nrow(crops)){
  if(is.list(ecocropPars(crops$NAME[i]))){
    crop <- ecocropPars(crops$NAME[i])
  } else {
    crop <- ecocropPars(crops$SCIENTNAME[i])
  }
  test2[i] <- crop$name %in% ecodb$nms
}
sink()
all(test2)
rm(test, test2, crop, i)

# RUN model
for(i in 1:nrow(crops)){
  if(is.list(ecocropPars(crops$NAME[i]))){
    crop <- ecocropPars(crops$NAME[i])
  } else {
    crop <- ecocropPars(crops$SCIENTNAME[i])
  }
  # add annual precipitation parameter
  ri <- which(crop$name == ecodb$nms)
  np <- cbind(annprec = as.numeric(ecodb[ri, raincols]))
  crop$parameters <- cbind(crop$parameters, np)
  
  m <- ecocrop(crop)
  control(m, get_max=TRUE)
  fn <- paste0("ID", formatC(crops$RID[i], width = 4, flag = "0"))
  predict(m, tavg = tavg, ktmp = tmin, prec = prec, ph = ph, annprec = annprec,
          filename = file.path(dirn, paste0(fn, ".tif")), 
          overwrite = T, 
          wopt = list(names = fn, filetype = "GTiff", 
                      gdal = c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))
}

