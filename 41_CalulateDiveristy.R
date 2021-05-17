library(magrittr)
library(data.table)
library(terra)
library(rnaturalearth)
library(viridis)

if(!grepl("globcropdiv$", getwd())){
  if(system('hostname', TRUE) %in% c("ESP-RH-9891", "LAPTOP-ST129J47")){ 
    setwd("G:/My Drive/globcropdiv/")
  } # else if { ... 
}

# diversity function 
fd <- function(x, na.rm = na.rm){
  if(na.rm == TRUE){
    x <- x[x > 0 & !is.na(x)]  
  } else {
    x <- x[x > 0]
  }
  return(exp(-sum(x * log(x))))
}


# total cropland and countries
totcl <- rast("InData/TotalCropland.tif")

#countries <- ne_download(scale = 10, type = "countries")
countries <- vect("InData/countries/ne_10m_admin_0_countries.shp")


# Actual ########################################
actual <- "InData/CropAbundance/*.tif" %>%  Sys.glob() %>%  rast()
crops <- names(actual)

d <- values(totcl) %>% 
  {data.table(cell = which(!is.na(.)), tcl = .[!is.na(.)])}
d[, (crops):= extract(actual, cell)]
d[, (crops):= lapply(.SD, `/`, tcl), .SDcols = crops]

d[, Da:= apply(.SD, 1, fd, na.rm = TRUE), .SDcols = crops]

# raster
r <- rast(totcl)
v <- rep(NA_real_, ncell(r))
v[d$cell] <- d$Da
values(r) <- v

writeRaster(r, filename = "OutData/Da.tif", overwrite = T, 
            wopt = list(names = "Da", filetype = "GTiff",
                        gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")                        )
)




# Potential ####################################
## Integrated Models ----

# Upload allocated Integrated data
d <- fread("OutData/allocated_int.csv") 
acols <- names(d)[names(d) %like% '^a\\.']

# total allocated area
d[, tcl:= Reduce(`+`, .SD), .SDcols = acols]

# allocated area fraction
fcols <- gsub("^a", "f", acols)
d[, (fcols):= .SD/tcl, .SDcols = acols]

# Potential Diversity
d[, Dp:= apply(.SD, 1, fd, na.rm = TRUE), .SDcols = fcols]

# raster
r <- rast(totcl)
v <- rep(NA_real_, ncell(r))
v[d$cell] <- d$Dp
values(r) <- v

writeRaster(r, filename = "OutData/Dp_int.tif", overwrite = T, 
            wopt = list(names = "Dp_int", filetype = "GTiff",
                        gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")                        )
)


## SDM Models ----

# Upload allocated SDM data
d <- fread("OutData/allocated_sdm.csv") 
acols <- names(d)[names(d) %like% '^a\\.']

# total allocated area
d[, tcl:= Reduce(`+`, .SD), .SDcols = acols]

# allocated area fraction
fcols <- gsub("^a", "f", acols)
d[, (fcols):= .SD/tcl, .SDcols = acols]

# Potential Diversity
d[, Dp:= apply(.SD, 1, fd, na.rm = TRUE), .SDcols = fcols]

# raster
r <- rast(totcl)
v <- rep(NA_real_, ncell(r))
v[d$cell] <- d$Dp
values(r) <- v

writeRaster(r, filename = "OutData/Dp_sdm.tif", overwrite = T, 
            wopt = list(names = "Dp_sdm", filetype = "GTiff",
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

# Potential Diversity
d[, Dp:= apply(.SD, 1, fd, na.rm = TRUE), .SDcols = fcols]

# raster
r <- rast(totcl)
v <- rep(NA_real_, ncell(r))
v[d$cell] <- d$Dp
values(r) <- v

writeRaster(r, filename = "OutData/Dp_eco.tif", overwrite = T, 
            wopt = list(names = "Dp_eco", filetype = "GTiff",
                        gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")                        )
)

