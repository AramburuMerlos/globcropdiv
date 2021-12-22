library(terra)

gdal()
#### THIS SCRIPT NEEDS GDAL >=3.2.2 !!!!!!!!!!


library(data.table)

setwd("D:/globcropdiv/")

dir <- "OutData/projected"
dir.create(dir, F, T)

# model raster
totcl <- rast("InData/TotalCropland.tif")
r <- rast(totcl)
rproj <- project(r, "+proj=eqearth +ellps=WGS84 +datum=WGS84")
res(rproj) <- 9260


# Total Cropland ##########################
# project proportions instead of areas
cell_area <- cellSize(totcl)/10000 # ha
cl_prop <- totcl/cell_area

project(cl_prop, rproj, 
        filename = file.path(dir, "CroplandProp.tif"), overwrite = T, 
        wopt = list(names = "tcl_prop", filetype = "GTiff",
                    gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))


# Actual ########################################

actual <- rast(Sys.glob("InData/CropAbundance/*.tif"))
actual_prop <- actual/totcl

project(actual_prop, rproj, 
        filename = file.path(dir, "ActualCropProp.tif"), overwrite = T, 
        wopt = list(names = names(actual), filetype = "GTiff",
                    gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))



# Potential ##################################

## Ecocrop #####################################

suit <- rast(Sys.glob("OutData/Ecocrop/*.tif"))
totsuit <- app(suit, sum, na.rm = T)
pot_eco_prop <- suit/totsuit

project(pot_eco_prop, rproj, 
        filename = file.path(dir, "PotEcoCropProp.tif"), overwrite = T, 
        wopt = list(names = names(suit), filetype = "GTiff",
                    gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))


## SDM #####################################
suit <- rast(Sys.glob("OutData/SDM/*AVG.tif"))
crops <- names(suit)

# because of memory issues
tdir <- "TmpData/pot_sdm_prop"
dir.create(tdir, F, T)

app(suit, sum, na.rm = T, 
    filename = "TmpData/SumOfSDMSuit.tif", overwrite = T, 
    wopt = list(names = "totsuit", filetype = "GTiff",
                gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))
gc(reset = T)

totsuit <- rast("TmpData/SumOfSDMSuit.tif")
vts <- values(totsuit)[,1]

for(i in 1:length(crops)){
  r <- rast(totsuit)
  v <- values(suit[[i]])[,1]
  v <- v/vts
  values(r) <- v
  
  writeRaster(
    r, overwrite = T, 
    filename = file.path(tdir, paste0(crops[i], ".tif")), 
    wopt = list(names = crops[i], filetype = "GTiff",
                gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")                        )
  )
}

pot_sdm_prop <- rast(Sys.glob(paste0(tdir, "/*.tif")))

project(pot_sdm_prop, rproj, 
        filename = file.path(dir, "PotSDMCropProp.tif"), overwrite = T, 
        wopt = list(names = names(suit), filetype = "GTiff",
                    gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))




# Attainable ####################################

## SDM ----

# Upload allocated SDM data
d <- fread("OutData/allocated_sdm.csv") 
acols <- names(d)[names(d) %like% '^a\\.']
crops <- gsub("^a\\.", "", acols)

d[, tcl:= Reduce(`+`, .SD), .SDcols = acols]
d[, (acols):= lapply(.SD, function(x) x/tcl), .SDcols = acols]

tdir <- "TmpData/att_sdm_prop"
dir.create(tdir, F, T)

for(i in 1:length(crops)){
  r <- rast(totcl)
  v <- rep(NA_real_, ncell(r))
  v[d$cell] <- d[[acols[i]]]
  values(r) <- v
  
  writeRaster(
    r, overwrite = T, 
    filename = file.path(tdir, paste0(crops[i], ".tif")), 
    wopt = list(names = crops[i], filetype = "GTiff",
                gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")                        )
  )
}

rm(d)

att_sdm_prop <- rast(Sys.glob(paste0(tdir, "/*.tif")))

project(att_sdm_prop, rproj, 
        filename = file.path(dir, "AttSDMCropProp.tif"), overwrite = T, 
        wopt = list(names = names(att_sdm_prop), filetype = "GTiff",
                    gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))



## Ecocrop ----

# Upload allocated SDM data
d <- fread("OutData/allocated_eco.csv") 
acols <- names(d)[names(d) %like% '^a\\.']
crops <- gsub("^a\\.", "", acols)

d[, tcl:= Reduce(`+`, .SD), .SDcols = acols]
d[, (acols):= lapply(.SD, function(x) x/tcl), .SDcols = acols]

tdir <- "TmpData/att_eco_prop"
dir.create(tdir, F, T)

for(i in 1:length(crops)){
  r <- rast(totcl)
  v <- rep(NA_real_, ncell(r))
  v[d$cell] <- d[[acols[i]]]
  values(r) <- v
  
  writeRaster(
    r, overwrite = T, 
    filename = file.path(tdir, paste0(crops[i], ".tif")), 
    wopt = list(names = crops[i], filetype = "GTiff",
                gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")                        )
  )
}

rm(d)

att_eco_prop <- rast(Sys.glob(paste0(tdir, "/*.tif")))

project(att_eco_prop, rproj, 
        filename = file.path(dir, "AttEcoCropProp.tif"), overwrite = T, 
        wopt = list(names = names(att_eco_prop), filetype = "GTiff",
                    gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))


 
 
 # Sum of proj Crop prop #####
actual <- rast("OutData/projected/ActualCropProp.tif")
app(actual, sum, na.rm = TRUE) 

pot_eco <- rast("OutData/projected/PotEcoCropProp.tif")
app(pot_eco, sum, na.rm = TRUE) 

pot_sdm <- rast("OutData/projected/PotSDMCropProp.tif")
app(pot_sdm, sum, na.rm = TRUE) 

att_eco <- rast("OutData/projected/AttEcoCropProp.tif")
app(att_eco, sum, na.rm = TRUE) 

att_sdm <- rast("OutData/projected/AttSDMCropProp.tif")
app(att_sdm, sum, na.rm = TRUE) 


