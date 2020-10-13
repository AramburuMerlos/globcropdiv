library(terra)
library(Recocrop)

nwd <- nchar(getwd())
if(substr(getwd(),  nwd - 10, nwd) != "globcropdiv") stop("See 0000_wd.R")

if(all(sapply(c("D:/WorldClim", "D:/SoilGrids", "D:/globcropdiv"), dir.exists))){
  wcpath <- "D:/WorldClim"
  sgpath <- "D:/SoilGrids"
  gdpath <- "D:/globcropdiv"
} else {
  wcpath <- "InData/WorldClim"
  sgpath <- "InData/SoilGrids"
  gdpath <- "OutData"
}

# climatic data
tavg <- rast(Sys.glob(file.path(wcpath, "2.1/wc5min/tavg/*.tif")))
tmin <- rast(Sys.glob(file.path(wcpath, "2.1/wc5min/tmin/*.tif")))

# ph
ph <- rast(file.path(sgpath, "phh2o/phh2o_0-15cm_mean_5min.tif"))
ph <- ph/10

# crop mask
cm <- rast(file.path(gdpath, "Monf_Cropland_Mask.tif"))

# mask data
tavg <- mask(tavg, cm)
tmin <- mask(tmin, cm)
ph <- mask(ph, cm)

# crops to be modelled
dcrops <- read.csv(here::here("AuxData/Monfcrops.csv"), 
                  stringsAsFactors = F) 
crops <- unique(dcrops[, c("RID", "NAME", "SCIENTNAME")])
crops <- crops[order(crops$RID),]

dirn <- (paste0(gdpath, "/EcocropSuit/byEcoID/Irrigated"))
dir.create(dirn)
dirn <- (paste0(gdpath, "/EcocropSuit/byEcoID/Irrigated/RID"))
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
  predict(m, tavg = tavg, ktmp = tmin, ph = ph, 
          filename = file.path(dirn, paste0(fn, ".tif")), 
          overwrite = T, 
          wopt = list(names = fn, filetype = "GTiff",
          gdal = c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))
}

