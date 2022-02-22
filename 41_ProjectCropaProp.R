library(terra)

gdal()
#### THIS SCRIPT NEEDS GDAL >=3.2.2 !!!!!!!!!!


library(data.table)

if(system('hostname', TRUE) %in% c("ESP-RH-9891", "LAPTOP-ST129J47")){
  setwd("D:/globcropdiv/")
} # else if { ... 


dir <- "OutData/projected"
dir.create(dir, F, T)

# model raster
totcl <- rast("InData/TotalCropland.tif")
r <- rast(totcl)
rproj <- project(r, "+proj=eqearth",  mask = TRUE)
res(rproj) <- 9260


# Total Cropland ##########################
# project proportions instead of areas
cell_area <- cellSize(totcl)/10000 # ha
cl_prop <- totcl/cell_area

project(cl_prop, rproj,  mask = TRUE,
        filename = file.path(dir, "CroplandProp.tif"), overwrite = T, 
        wopt = list(names = "tcl_prop", filetype = "GTiff",
                    gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))


# Actual ########################################

fn <- Sys.glob("InData/CropAbundance/*.tif")
crops <- gsub("InData/CropAbundance/", "", gsub(".tif$", "", fn))

tdir <- "TmpData/actual_prop"
dir.create(tdir, F, T)

for(i in 1:length(crops)){
  rtmp <- rast(fn[i])
  rtmp <- rtmp/totcl
  
  writeRaster(
    rtmp, overwrite = T, 
    filename = file.path(tdir, paste0(crops[i], ".tif")), 
    wopt = list(names = crops[i], filetype = "GTiff",
                gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")                        )
  )
}

actual_prop <- rast(Sys.glob(paste0(tdir, "/*.tif")))

project(actual_prop, rproj,  mask = TRUE,
        filename = file.path(dir, "ActualCropProp.tif"), overwrite = T, 
        wopt = list(names = crops, filetype = "GTiff",
                    gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))



# Potential ##################################

## Ecocrop #####################################
suit <- rast(Sys.glob("OutData/Ecocrop/*.tif"))
crops <- names(suit)

# because of memory issues
tdir <- "TmpData/pot_eco_prop"
dir.create(tdir, F, T)

app(suit, sum, na.rm = T, 
    filename = "TmpData/SumOfEcoSuit.tif", overwrite = T, 
    wopt = list(names = "totsuit", filetype = "GTiff",
                gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))
gc(reset = T)

totsuit <- rast("TmpData/SumOfEcoSuit.tif")
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
                gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6"))
  )
}

pot_eco_prop <- rast(Sys.glob(paste0(tdir, "/*.tif")))

project(pot_eco_prop, rproj,  mask = TRUE,
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

project(pot_sdm_prop, rproj,  mask = TRUE,
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

project(att_sdm_prop, rproj, mask = TRUE,
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

project(att_eco_prop, rproj, mask = TRUE, 
        filename = file.path(dir, "AttEcoCropProp.tif"), overwrite = T, 
        wopt = list(names = names(att_eco_prop), filetype = "GTiff",
                    gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))


 
 
# Sum of proj Crop prop ######
sdir <- "OutData/projected/sums"
dir.create(sdir, F, T)

actual <- rast("OutData/projected/ActualCropProp.tif")

app(
  actual, sum, na.rm = TRUE, 
  filename = file.path(sdir, "SumActualCropProp.tif"), overwrite = T, 
  wopt = list(filetype = "GTiff", 
              gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6"))
)
 

pot_eco <- rast("OutData/projected/PotEcoCropProp.tif")

app(
  pot_eco, sum, na.rm = TRUE,
  filename = file.path(sdir, "SumPotEcoCropProp.tif"), overwrite = T, 
  wopt = list(filetype = "GTiff",
              gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6"))
) 


pot_sdm <- rast("OutData/projected/PotSDMCropProp.tif")

app(
  pot_sdm, sum, na.rm = TRUE,
  filename = file.path(sdir, "SumPotSDMCropProp.tif"), overwrite = T, 
  wopt = list(filetype = "GTiff",
              gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6"))

) 


att_eco <- rast("OutData/projected/AttEcoCropProp.tif")

app(
  att_eco, sum, na.rm = TRUE,
  filename = file.path(sdir, "SumAttEcoCropProp.tif"), overwrite = T, 
  wopt = list(filetype = "GTiff",
              gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6"))
) 


att_sdm <- rast("OutData/projected/AttSDMCropProp.tif")

app(
  att_sdm, sum, na.rm = TRUE,
  filename = file.path(sdir, "SumAttSDMCropProp.tif"), overwrite = T, 
  wopt = list(filetype = "GTiff",
              gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6"))
) 


# Adjust prop to sum 1 ----------
dir <- "OutData/projected"
sdir <- "OutData/projected/sums"
dir.create(sdir, F, T)

actual <- rast("OutData/projected/ActualCropProp.tif")
sumactual <- rast(file.path(sdir, "SumActualCropProp.tif"))
sumactual <- classify(sumactual, cbind(0,NA))

tmp <- actual/sumactual

writeRaster(
  tmp,
  filename = file.path(dir, "ActualCropProp.tif"), overwrite = T, 
  wopt = list(names = names(actual), filetype = "GTiff",
              gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6"))
)
rm(tmp, actual, sumactual)
gc()


pot_eco <- rast("OutData/projected/PotEcoCropProp.tif")
sum_pot_eco <- rast(file.path(sdir, "SumPotEcoCropProp.tif"))
sum_pot_eco <- classify(sum_pot_eco, cbind(0, NA))

tmp <- pot_eco/sum_pot_eco

writeRaster(
  tmp,
  filename = file.path(dir, "PotEcoCropProp.tif"), overwrite = T, 
  wopt = list(names = names(pot_eco), filetype = "GTiff",
              gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6"))
)

rm(tmp, pot_eco, sum_pot_eco)
gc()


pot_sdm <- rast("OutData/projected/PotSDMCropProp.tif")
sum_pot_sdm <- rast(file.path(sdir, "SumPotSDMCropProp.tif"))
sum_pot_sdm <- classify(sum_pot_sdm, cbind(0, NA))

tmp <- pot_sdm/sum_pot_sdm

writeRaster(
  tmp,
  filename = file.path(dir, "PotSDMCropProp.tif"), overwrite = T, 
  wopt = list(names = names(pot_sdm), filetype = "GTiff",
              gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6"))
)
rm(tmp, pot_sdm, sum_pot_sdm)
gc()

att_eco <- rast("OutData/projected/AttEcoCropProp.tif")
sum_att_eco <- rast(file.path(sdir, "SumAttEcoCropProp.tif"))
sum_att_eco <- classify(sum_att_eco, cbind(0, NA))

tmp <- att_eco/sum_att_eco

writeRaster(
  tmp,
  filename = file.path(dir, "AttEcoCropProp.tif"), overwrite = T, 
  wopt = list(names = names(att_eco), filetype = "GTiff",
              gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6"))
)
rm(tmp, att_eco, sum_att_eco)
gc()

att_sdm <- rast("OutData/projected/AttSDMCropProp.tif")
sum_att_sdm <- rast(file.path(sdir, "SumAttSDMCropProp.tif"))
sum_att_sdm <- classify(sum_att_sdm, cbind(0, NA))

tmp <- att_sdm/sum_att_sdm

writeRaster(
  tmp,
  filename = file.path(dir, "AttSDMCropProp.tif"), overwrite = T, 
  wopt = list(names = names(att_sdm), filetype = "GTiff",
              gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6"))
)
rm(tmp, att_sdm, sum_att_sdm)
gc()



# crop proportion gap SDM -----------
dir <- "OutData/projected"

# crop proportions in low diversity cells
actual_area <- rast("OutData/projected/ActualCropProp.tif")
crops <- names(actual_area)

attain_area <-rast("OutData/projected/AttSDMCropProp.tif")
all.equal(names(attain_area), crops)

gap <- attain_area - actual_area

writeRaster(
  gap,
  filename = file.path(dir, "CropPropGap_SDM.tif"), overwrite = T, 
  wopt = list(names = names(attain_area), filetype = "GTiff",
              gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6"))
)

# crop proportion gap Eco -----------
dir <- "OutData/projected"

# crop proportions in low diversity cells
actual_area <- rast("OutData/projected/ActualCropProp.tif")
crops <- names(actual_area)

attain_area <-rast("OutData/projected/AttEcoCropProp.tif")
all.equal(names(attain_area), crops)

gap <- attain_area - actual_area

writeRaster(
  gap,
  filename = file.path(dir, "CropPropGap_Eco.tif"), overwrite = T, 
  wopt = list(names = names(attain_area), filetype = "GTiff",
              gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6"))
)



