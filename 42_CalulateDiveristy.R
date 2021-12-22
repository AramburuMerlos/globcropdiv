library(data.table)
library(terra)


if(system('hostname', TRUE) == "ESP-RH-9891"){
  setwd("D:/globcropdiv/")
} else if(system('hostname', TRUE) == "LAPTOP-ST129J47"){ 
  setwd("G:/My Drive/globcropdiv/")
} # else if { ... 

# diversity function #############
fd <- function(x, na.rm = na.rm){
  if(na.rm == TRUE){
    x <- x[x > 0 & !is.na(x)]  
  } else {
    x <- x[x > 0]
  }
  return(exp(-sum(x * log(x))))
}

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

# raster
r <- rast(cl_prop)

# Actual D ########################################
d <- data.table(cell = which(!is.na(values(cl_prop))))

d[, (crops):= extract(actual, cell)]
d[, D:= apply(.SD, 1, fd, na.rm = TRUE), .SDcols = crops]

v <- rep(NA_real_, ncell(r))
v[d$cell] <- d$D
values(r) <- v
writeRaster(r, 
            filename = "OutData/act_D.tif", overwrite = T, 
            wopt = list(names = "act_D", filetype = "GTiff",
                        gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))


# Potential ##################################

## Ecocrop #####################################
d[, (crops):= extract(pot_eco, cell)]
d[, D:= apply(.SD, 1, fd, na.rm = TRUE), .SDcols = crops]

# raster
v <- rep(NA_real_, ncell(r))
v[d$cell] <- d$D
values(r) <- v

writeRaster(r, filename = "OutData/pot_D_eco.tif", overwrite = T, 
            wopt = list(names = "pot_D_eco", filetype = "GTiff",
                        gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))

## SDM #####################################
d[, (crops):= extract(pot_sdm, cell)]
d[, D:= apply(.SD, 1, fd, na.rm = TRUE), .SDcols = crops]

v <- rep(NA_real_, ncell(r))
v[d$cell] <- d$D
values(r) <- v

writeRaster(r, filename = "OutData/pot_D_sdm.tif", overwrite = T, 
            wopt = list(names = "pot_D_sdm", filetype = "GTiff",
                        gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))


# Attainable ##################################

## Ecocrop #####################################
d[, (crops):= extract(att_eco, cell)]
d[, D:= apply(.SD, 1, fd, na.rm = TRUE), .SDcols = crops]

# raster
v <- rep(NA_real_, ncell(r))
v[d$cell] <- d$D
values(r) <- v

writeRaster(r, filename = "OutData/att_D_eco.tif", overwrite = T, 
            wopt = list(names = "att_D_eco", filetype = "GTiff",
                        gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))

## SDM #####################################
d[, (crops):= extract(att_sdm, cell)]
d[, D:= apply(.SD, 1, fd, na.rm = TRUE), .SDcols = crops]

# raster
v <- rep(NA_real_, ncell(r))
v[d$cell] <- d$D
values(r) <- v

writeRaster(r, filename = "OutData/att_D_sdm.tif", overwrite = T, 
            wopt = list(names = "att_D_sdm", filetype = "GTiff",
                        gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))
