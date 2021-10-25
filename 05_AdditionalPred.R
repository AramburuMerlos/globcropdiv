library(terra)
library(geosphere)
library(data.table)

if(system('hostname', TRUE) == "ESP-RH-9891"){
  setwd("D:/globcropdiv/")
} else if(system('hostname', TRUE) == "LAPTOP-ST129J47"){ 
  setwd("G:/My Drive/globcropdiv/")
} # else if { ... 

r <- rast("InData/TotalCropland.tif")

# Elevation heterogeneity -----------
# GMTED2010 downloaded from 
# https://edcintl.cr.usgs.gov/downloads/sciweb1/shared/topo/downloads/GMTED/
r_elev_sd <- rast("InData/GMTED2010/sd30_grd/w001001.adf")

# fix geometry
r_elev_sd <- extend(r_elev_sd, r)
f_agg_sd <- function(x) sqrt(mean(x^2, na.rm = T))
r_elev_sd <- aggregate(r_elev_sd, fact = 10, fun = f_agg_sd)

compareGeom(r, r_elev_sd)

writeRaster(r_elev_sd, overwrite = T, 
            filename = "InData/GMTED2010/elev_sd.tif", 
            wopt = list(names = "field_size", filetype = "GTiff",
                        gdal = c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))


# Field size ---------
# Data from Friz et al., 2015 (https://doi.org/10.1111/gcb.12838)
# downloaded from geo wiki app
r_field <- rast(Sys.glob("InData/field_size*/field_size*.img"))


## Fix geometry ------------
# change 0 field size to NA so it doesn't affect the aggregation
r_field <- extend(r_field, r)
r_field <- classify(r_field, cbind(0, NA))
r_field <- aggregate(r_field, fact = 10, fun = "median", na.rm = TRUE)

compareGeom(r, r_field)


## Fill NA -------------------

# field data has NA values from cells where there is cropland. 
# data will be filled by idw 
f_idw <- function(x, k = 2){
  # one dimension distance to center
  xc <- ceiling(x/2)
  d1 <- abs(1:x - xc)
  md1 <- matrix(rep(d1, x), nrow = x, ncol = x)
  md2 <- t(md1)
  m <- sqrt(md1^2 + md2^2)
  idw <- (1/(m^k))
  idw[is.infinite(idw)] <- NA
  idw <- idw/sum(idw, na.rm = T)
  return(idw)
}


d <- data.table(cell = which(!is.na(values(r))))
d[, field_size:= extract(r_field, cell)]

n_na <- sum(is.na(d$field_size))
x <- 1

while(n_na > 0){    
  x <- x + 10
  w <- f_idw(x)
  
  r_field <- focal(r_field, w = w, fun = "mean", na.only = TRUE)
  # the focal returns 0 when is all NA. 
  r_field <- classify(r_field, cbind(0, NA))
  r_field <- mask(r_field, r) # it fills everything, including oceans
  d[is.na(field_size), field_size:= extract(r_field, cell)]
  n_na <- sum(is.na(d$field_size))
  print(paste("matrix side:", x, "; remaining NA:", n_na))
  if(x > 50) break # the focal function gets stuck with a big matrix
}

# let's try a different approach for the very isolated NA's

w_na <- which(is.na(d$field_size))

d[, `:=`(lon = xFromCell(r, cell), lat = yFromCell(r, cell))]

# distance from i cell to cells less than x degrees apart (N-S or E-W)
x <- 10
while(length(w_na) > 0){
  for(i in w_na){
    scells <- which(abs(d[i, lon] - d[,lon]) < x, abs(d[i, lat] - d[,lat]) < x)  
    d[scells, gdist:= distGeo(d[i, c(lon,lat)], cbind(lon, lat))]
    d[scells, w:= 1/(gdist^2)]
    d[scells, w:= w/sum(w, na.rm = T)]
    d[i, field_size:= d[scells, sum(field_size * w, na.rm = T)]]
  }
  w_na <- which(is.na(d$field_size))
  x <- x + 10
}

v <- rep(NA_real_, ncell(r_field))
v[d$cell] <- d$field_size

values(r_field) <- v

writeRaster(r_field, overwrite = T, 
            filename = "InData/field_size_10_40_cropland/field_size.tif", 
            wopt = list(names = "field_size", filetype = "GTiff",
                        gdal = c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))







