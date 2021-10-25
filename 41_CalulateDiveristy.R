library(magrittr)
library(data.table)
library(terra)


if(system('hostname', TRUE) == "ESP-RH-9891"){
  setwd("D:/globcropdiv/")
} else if(system('hostname', TRUE) == "LAPTOP-ST129J47"){ 
  setwd("G:/My Drive/globcropdiv/")
} # else if { ... 

# diversity function 
fd <- function(x, na.rm = na.rm){
  if(na.rm == TRUE){
    x <- x[x > 0 & !is.na(x)]  
  } else {
    x <- x[x > 0]
  }
  return(exp(-sum(x * log(x))))
}


# total cropland 
totcl <- rast("InData/TotalCropland.tif")


# Actual ########################################

actual <- "InData/CropAbundance/*.tif" %>%  Sys.glob() %>%  rast()
crops <- names(actual)

d <- values(totcl) %>% 
  {data.table(cell = which(!is.na(.)), tcl = .[!is.na(.)])}
d[, (crops):= extract(actual, cell)]
d[, (crops):= lapply(.SD, `/`, tcl), .SDcols = crops]

d[, act_D:= apply(.SD, 1, fd, na.rm = TRUE), .SDcols = crops]

# raster
r <- rast(totcl)
v <- rep(NA_real_, ncell(r))
v[d$cell] <- d$act_D
values(r) <- v

writeRaster(r, filename = "OutData/act_D.tif", overwrite = T, 
            wopt = list(names = "act_D", filetype = "GTiff",
                        gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")                        )
)


# Potential ##################################

## Ecocrop #####################################

suit <- "OutData/Ecocrop/*.tif" %>%  Sys.glob() %>%  rast()
crops <- names(suit)

d <- data.table(cell = which(!is.na(values(totcl))))
d[, (crops):= extract(suit, cell)]

d[, rs:= Reduce(`+`, .SD), .SDcols = crops]
d[, (crops):= lapply(.SD, `/`, rs), .SDcols = crops]

d[, D:= apply(.SD, 1, fd, na.rm = TRUE), .SDcols = crops]

# raster
r <- rast(totcl)
v <- rep(NA_real_, ncell(r))
v[d$cell] <- d$D
values(r) <- v

writeRaster(r, filename = "OutData/pot_D_eco.tif", overwrite = T, 
            wopt = list(names = "pot_D_eco", filetype = "GTiff",
                        gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")                        )
)

## SDM #####################################
suit <- "OutData/SDM/*AVG.tif" %>%  Sys.glob() %>%  rast()
crops <- names(suit)

d <- data.table(cell = which(!is.na(values(totcl))))
d[, (crops):= extract(suit, cell)]

d[, rs:= Reduce(`+`, .SD), .SDcols = crops]
d[, (crops):= lapply(.SD, `/`, rs), .SDcols = crops]

d[, D:= apply(.SD, 1, fd, na.rm = TRUE), .SDcols = crops]

# raster
r <- rast(totcl)
v <- rep(NA_real_, ncell(r))
v[d$cell] <- d$D
values(r) <- v

writeRaster(r, filename = "OutData/pot_D_sdm.tif", overwrite = T, 
            wopt = list(names = "pot_D_sdm", filetype = "GTiff",
                        gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")                        )
)


# Attainable ####################################
## SDM Models ----

# Upload allocated SDM data
d <- fread("OutData/allocated_sdm.csv") 
acols <- names(d)[names(d) %like% '^a\\.']

# total allocated area
d[, tcl:= Reduce(`+`, .SD), .SDcols = acols]

# allocated area fraction
fcols <- gsub("^a", "f", acols)
d[, (fcols):= .SD/tcl, .SDcols = acols]

# Attainable Diversity
d[, D:= apply(.SD, 1, fd, na.rm = TRUE), .SDcols = fcols]

# raster
r <- rast(totcl)
v <- rep(NA_real_, ncell(r))
v[d$cell] <- d$D
values(r) <- v

writeRaster(r, filename = "OutData/att_D_sdm.tif", overwrite = T, 
            wopt = list(names = "att_D_sdm", filetype = "GTiff",
                        gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")                        )
)


## Ecocrop Model ----

# Upload allocated Ecocrop data
d <- fread("OutData/allocated_eco.csv") 
acols <- names(d)[names(d) %like% '^a\\.']

# total allocated area
d[, tcl:= Reduce(`+`, .SD), .SDcols = acols]

# allocated area fraction
fcols <- gsub("^a", "f", acols)
d[, (fcols):= .SD/tcl, .SDcols = acols]

# Attainable Diversity
d[, D:= apply(.SD, 1, fd, na.rm = TRUE), .SDcols = fcols]

# raster
r <- rast(totcl)
v <- rep(NA_real_, ncell(r))
v[d$cell] <- d$D
values(r) <- v

writeRaster(r, filename = "OutData/att_D_eco.tif", overwrite = T, 
            wopt = list(names = "att_D_eco", filetype = "GTiff",
                        gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")                        )
)

