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

# Data prep ---------

# Upload allocated SDM data
d <- fread("OutData/allocated_sdm.csv") 

crops <- gsub("^a\\.", "", names(d)[names(d) %like% '^a\\.'])
crops_allo_ha <- paste0(crops, "_allo_ha")
crops_allo_fr <- paste0(crops, "_allo_fr")
crops_actu_ha <- paste0(crops, "_actu_ha")
crops_actu_fr <- paste0(crops, "_actu_fr")
crops_suit <- paste0(crops, "_suit")

setnames(d, paste0("a.", crops), crops_allo_ha)
setnames(d, paste0("s.", crops), crops_suit)

# Actual Area
"InData/CropAbundance/*.tif" %>% 
  Sys.glob() %>% 
  rast() %>% 
  {d[, (crops_actu_ha):= extract(., cell, list = TRUE)]}

# actual area fraction
d[, (crops_actu_fr):= .SD/tcl, .SDcols = crops_actu_ha]

# total allocated area
d[, tacl:= Reduce(`+`, .SD), .SDcols = crops_allo_ha]

# allocated area fraction
d[, (crops_allo_fr):= .SD/tacl, .SDcols = crops_allo_ha]

# Calculate Diversity ------------
 
# diversity function
f <- function(x){
  x <- x[x > 0]
  return(exp(-sum(x * log(x))))
}

# Actual Diversity 
d[, Da:= apply(.SD, 1, f), .SDcols = crops_actu_fr]

# Potential Diversity
d[, Dp:= apply(.SD, 1, f), .SDcols = crops_allo_fr]

# Mapping  ###################################
countries <- ne_download(scale = 10, type = "countries")

r_tcl <- rast("InData/TotalCropland.tif")
r <- rast(r_tcl)
v <- rep(NA_real_, ncell(r))

ar = ncol(r)/nrow(r)
dir.create("Maps", F, T)

## Diversity -------------

### Actual ------
v[d$cell] <- log(d$Da)
values(r) <- v

tiff(filename = 'Maps/Da.tif', units = "in",
     width = ncol(r)/300, 
     height = (ncol(r)/300)/ar, 
     type = "cairo", res = 300, 
     compression = "zip")
par(mai = c(0,0,0,0), omi = c(0,0,0,0))
plot(r, type = "continuous", axes = FALSE, 
     col = turbo(256, direction = -1),
     maxcell = ncell(r), 
     mar = c(2,0,2,4), main = "Actual Diversity Index (Shannon)")
plot(countries, lwd = 0.4, add = T)
dev.off()

### Potential ------
v[d$cell] <- log(d$Dp)
values(r) <- v

tiff(filename = 'Maps/Dp_sdm.tif', units = "in",
     width = ncol(r)/300, 
     height = (ncol(r)/300)/ar, 
     type = "cairo", res = 300, 
     compression = "zip")
par(mai = c(0,0,0,0), omi = c(0,0,0,0))
plot(r, type = "continuous", axes = FALSE, 
     col = turbo(256, direction = -1),
     maxcell = ncell(r), 
     mar = c(2,0,2,4), main = "Potential Diversity Index (Shannon)")
plot(countries, lwd = 0.4, add = T)
dev.off()


### Gap ----
v[d$cell] <- pmax(1 - (d$Da/d$Dp), 0)
values(r) <- v

tiff(filename = 'Maps/Dg.tif', units = "in",
     width = ncol(r)/300, 
     height = (ncol(r)/300)/ar, 
     type = "cairo", res = 300, 
     compression = "zip")
par(mai = c(0,0,0,0), omi = c(0,0,0,0))
plot(r, type = "continuous", axes = FALSE, 
     col = turbo(256),
     maxcell = ncell(r), 
     mar = c(2,0,2,4), main = "Diversity gap")
plot(countries, lwd = 0.4, add = T)
dev.off()


## Crop Areas ------------

### Actual -----
dir.create("Maps/ActualCropArea", F, T)

for(i in 1:length(crops)){
  v[d$cell] <- d[[crops_actu_fr[i]]]
  v[v == 0] <- NA_real_
  values(r) <- v

  tiff(filename = paste0('Maps/ActualCropArea/', crops[i], ".tif"), 
       units = "in",
       width = ncol(r)/300, 
       height = (ncol(r)/300)/ar, 
       type = "cairo", res = 300, 
       compression = "zip")
  par(mai = c(0,0,0,0), omi = c(0,0,0,0))
  plot(r, type = "continuous", axes = FALSE,
       col = mako(256, direction = -1),
       maxcell = ncell(r), mar = c(2,0,2,6), 
       main = paste(crops[i], "actual area (fraction)"))
  plot(countries, lwd = 0.4, add = T)
  dev.off()
}

### Allocated -----
dir.create("Maps/AllocatedCropArea/SDM", F, T)

for(i in 1:length(crops)){
  v[d$cell] <- d[[crops_allo_fr[i]]]
  v[v == 0] <- NA_real_
  values(r) <- v
  
  tiff(filename = paste0('Maps/AllocatedCropArea/SDM/', crops[i], ".tif"), 
       units = "in",
       width = ncol(r)/300, 
       height = (ncol(r)/300)/ar, 
       type = "cairo", res = 300, 
       compression = "zip")
  par(mai = c(0,0,0,0), omi = c(0,0,0,0))
  plot(r, type = "continuous", axes = FALSE,
       col = mako(256, direction = -1),
       maxcell = ncell(r), mar = c(2,0,2,6), 
       main = paste(crops[i], "allocated area (fraction)"))
  plot(countries, lwd = 0.4, add = T)
  dev.off()
}


## Suitability -----------
### Integrated -----------

dir.create("Maps/Suitability/Integrated", F, T)

for(i in 1:length(crops)){
  v[d$cell] <- d[[crops_suit[i]]]
  v[v == 0] <- NA_real_
  values(r) <- v
  
  tiff(filename = paste0('Maps/Suitability/Integrated/', crops[i], ".tif"), 
       units = "in",
       width = ncol(r)/300, 
       height = (ncol(r)/300)/ar, 
       type = "cairo", res = 300, 
       compression = "zip")
  par(mai = c(0,0,0,0), omi = c(0,0,0,0))
  plot(r, type = "continuous", axes = FALSE,
       col = mako(256, direction = -1),
       
       maxcell = ncell(r), mar = c(2,0,2,6), 
       main = paste(crops[i], "suitability (Int)"))
  plot(countries, lwd = 0.4, add = T)
  dev.off()
}
