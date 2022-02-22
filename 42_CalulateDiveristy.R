library(terra)

if(system('hostname', TRUE) %in% c("ESP-RH-9891", "LAPTOP-ST129J47")){
  setwd("D:/globcropdiv/")
} # else if { ... 


# diversity functions #############
fd <- function(x) exp(-sum(x * log(x), na.rm = T))
fd2 <- function(x) 1/sum(x^2, na.rm = T)


# uplaod proportions #########
actual <- rast("OutData/projected/ActualCropProp.tif")
pot_eco <- rast("OutData/projected/PotEcoCropProp.tif")
pot_sdm <- rast("OutData/projected/PotSDMCropProp.tif")
att_eco <- rast("OutData/projected/AttEcoCropProp.tif")
att_sdm <- rast("OutData/projected/AttSDMCropProp.tif")

cl_prop <- rast("OutData/projected/CroplandProp.tif")

## check crops #####
crops <- names(actual)
all(
  all.equal(crops, names(pot_eco)),
  all.equal(crops, names(pot_sdm)),
  all.equal(crops, names(att_eco)),
  all.equal(crops, names(att_sdm))
)



# Calculate D ########################################
# and mask non cropland
act_D <- app(actual, fd)
mask(
  act_D, cl_prop,
  filename = "OutData/act_D.tif", overwrite = T, 
  wopt = list(names = "act_D", filetype = "GTiff",
              gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6"))
)

pot_D_eco <- app(pot_eco, fd)
mask(
  pot_D_eco, cl_prop,  
  filename = "OutData/pot_D_eco.tif", overwrite = T, 
  wopt = list(names = "pot_D_eco", filetype = "GTiff",
              gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6"))
)

pot_D_sdm <- app(pot_sdm, fd)
mask(
  pot_D_sdm, cl_prop,
  filename = "OutData/pot_D_sdm.tif", overwrite = T, 
  wopt = list(names = "pot_D_sdm", filetype = "GTiff",
              gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6"))
)

att_D_eco <- app(att_eco, fd)
mask(
  att_D_eco, cl_prop,
  filename = "OutData/att_D_eco.tif", overwrite = T, 
  wopt = list(names = "att_D_eco", filetype = "GTiff",
              gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6"))
)

att_D_sdm <- app(att_sdm, fd)
mask(
  att_D_sdm, cl_prop,
  filename = "OutData/att_D_sdm.tif", overwrite = T, 
  wopt = list(names = "att_D_sdm", filetype = "GTiff",
              gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6"))
)


# Calculate D2 ########################################
# and mask non cropland
act_D2 <- app(actual, fd2)
act_D2 <- classify(act_D2, cbind(Inf, NA))
mask(
  act_D2, cl_prop,
  filename = "OutData/act_D2.tif", overwrite = T, 
  wopt = list(names = "act_D2", filetype = "GTiff",
              gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6"))
)

pot_D2_eco <- app(pot_eco, fd2)
pot_D2_eco <- classify(pot_D2_eco, cbind(Inf, NA))
mask(
  pot_D2_eco, cl_prop,  
  filename = "OutData/pot_D2_eco.tif", overwrite = T, 
  wopt = list(names = "pot_D2_eco", filetype = "GTiff",
              gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6"))
)

pot_D2_sdm <- app(pot_sdm, fd2)
pot_D2_sdm <- classify(pot_D2_sdm, cbind(Inf, NA))
mask(
  pot_D2_sdm, cl_prop,
  filename = "OutData/pot_D2_sdm.tif", overwrite = T, 
  wopt = list(names = "pot_D2_sdm", filetype = "GTiff",
              gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6"))
)

att_D2_eco <- app(att_eco, fd2)
att_D2_eco <- classify(att_D2_eco, cbind(Inf, NA))
mask(
  att_D2_eco, cl_prop,
  filename = "OutData/att_D2_eco.tif", overwrite = T, 
  wopt = list(names = "att_D2_eco", filetype = "GTiff",
              gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6"))
)

att_D2_sdm <- app(att_sdm, fd2)
att_D2_sdm <- classify(att_D2_sdm, cbind(Inf, NA))
mask(
  att_D2_sdm, cl_prop,
  filename = "OutData/att_D2_sdm.tif", overwrite = T, 
  wopt = list(names = "att_D2_sdm", filetype = "GTiff",
              gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6"))
)

