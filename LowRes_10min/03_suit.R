library(here)
library(terra)
library(Recocrop)

# climatic data
tavg <- rast(Sys.glob("D:/WorldClim/2.1/wc10min/tavg/*.tif")) 
tmin <- rast(Sys.glob("D:/WorldClim/2.1/wc10min/tmin/*.tif")) 
prec <- rast(Sys.glob("D:/WorldClim/2.1/wc10min/prec/*.tif")) 

# ph
ph <- rast("D:/SoilGrids/phh2o/phh20_0-5cm_10min.tif")
ph <- ph/10

# crop mask
cm <- rast("D:/GlobalCropMask_NASA/CropMask10min.tif")

# mask data
tavg <- mask(tavg, cm, maskvalue = 0)
tmin <- mask(tmin, cm, maskvalue = 0)
prec <- mask(prec, cm, maskvalue = 0)
ph <- mask(ph, cm, maskvalue = 0)

# crops to be modelled
crops <- read.csv(here("AuxData/selsp.csv"), stringsAsFactors = F) 

dir1 <- ("D:/Dp_world")
dir2 <- ("LowRes_10min")
dir3 <- ("CropSuit")

dir.create(dir1)
dir.create(file.path(dir1,dir2))
dir.create(file.path(dir1,dir2, dir3))

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
          filename = file.path(dir1,dir2,dir3, paste0(fn, ".tif")), 
          overwrite = T, 
          wopt = list(names = fn, filetype = "GTiff"))
}

