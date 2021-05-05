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

dir_area <- 'Maps/Area'
dir.create(dir_area, F, T)

dir_suit <- 'Maps/Suit'
dir.create(dir_suit, F, T)

#countries <- ne_download(scale = 10, type = "countries")
countries <- vect("InData/countries/ne_10m_admin_0_countries.shp")

# Total Cropland ----------------
totcl <- rast("InData/TotalCropland.tif")
ar = ncol(totcl)/nrow(totcl)

tiff(filename = paste0('Maps/Total_Cropland.tif'), 
     units = "in",
     width = ncol(totcl)/300, 
     height = (ncol(totcl)/300)/ar, 
     type = "cairo", res = 300, 
     compression = "zip")
par(mai = c(0,0,0,0), omi = c(0,0,0,0))
plot(totcl, type = "continuous", axes = FALSE,
     col = mako(256, direction = -1),
     maxcell = ncell(totcl), mar = c(2,0,2,6), 
     main = "Cropland (ha)")
plot(countries, lwd = 0.4, add = T)
dev.off()


# Actual  ########################################

## crop area -----
actual <- "InData/CropAbundance/*.tif" %>%  Sys.glob() %>%  rast()
crops <- names(actual)

d <- values(totcl) %>% 
  {data.table(cell = which(!is.na(.)), tcl = .[!is.na(.)])}
d[, (crops):= extract(actual, cell)]
d[, (crops):= lapply(.SD, `/`, tcl), .SDcols = crops]

# set 0 crop area as NA for plotting
for(j in crops) set(d, which(d[[j]] == 0), j, NA)

# empty raster
r <- rast(totcl)
v <- rep(NA_real_, ncell(r))


for(i in 1:length(crops)){
  v[d$cell] <- d[[crops[i]]]
  values(r) <- v
  
  tiff(filename = file.path(dir_area, paste0(crops[i], "_act.tif")), 
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


## Diversity -------------
r <- rast("OutData/Da.tif")
r <- log(r)

tiff(filename = 'Maps/Da.tif', units = "in",
     width = ncol(r)/300, 
     height = (ncol(r)/300)/ar, 
     type = "cairo", res = 300, 
     compression = "zip")
par(mai = c(0,0,0,0), omi = c(0,0,0,0))
plot(r, type = "continuous", axes = FALSE, 
     col = turbo(256, direction = -1),
     maxcell = ncell(r), 
     mar = c(2,0,2,4), main = "Actual Diversity Index")
plot(countries, lwd = 0.4, add = T)
dev.off()



# Integrated Models ####################################

## Data prep ---------
# Upload allocated SDM data
d <- fread("OutData/allocated_int.csv") 

acols <- names(d)[names(d) %like% '^a\\.']
scols <- gsub("^a", "^s", acols)
crops <- gsub("^a\\.", "", acols)


# total allocated area
d[, tcl:= Reduce(`+`, .SD), .SDcols = acols]

# allocated area fraction
d[, (acols):= .SD/tcl, .SDcols = acols]

# set 0 crop area as NA for plotting
for(j in acols) set(d, which(d[[j]] == 0), j, NA)

# set 0 crop suit as NA for plotting
for(j in scols) set(d, which(d[[j]] == 0), j, NA)

# empty raster
r <- rast(totcl)
v <- rep(NA_real_, ncell(r))

## crop area ------------
for(i in 1:length(acols)){
  v[d$cell] <- d[[acols[i]]]
  values(r) <- v
  
  tiff(filename = file.path(dir_area, paste0(crops[i], "_int.tif")),  
       units = "in",
       width = ncol(r)/300, 
       height = (ncol(r)/300)/ar, 
       type = "cairo", res = 300, 
       compression = "zip")
  par(mai = c(0,0,0,0), omi = c(0,0,0,0))
  plot(r, type = "continuous", axes = FALSE,
       col = mako(256, direction = -1),
       maxcell = ncell(r), mar = c(2,0,2,6), 
       main = paste(crops[i], "allocated area (int)"))
  plot(countries, lwd = 0.4, add = T)
  dev.off()
}

## suitability -----------
for(i in 1:length(scols)){
  v[d$cell] <- d[[scols[i]]]
  values(r) <- v
  
  tiff(filename = file.path(dir_suit, paste0(crops[i], "_int.tif")),  
       units = "in",
       width = ncol(r)/300, 
       height = (ncol(r)/300)/ar, 
       type = "cairo", res = 300, 
       compression = "zip")
  par(mai = c(0,0,0,0), omi = c(0,0,0,0))
  plot(r, type = "continuous", axes = FALSE,
       col = mako(256, direction = -1),
       maxcell = ncell(r), mar = c(2,0,2,6), 
       main = paste(crops[i], "suitability (int)"))
  plot(countries, lwd = 0.4, add = T)
  dev.off()
}

## Diversity -------------
r <- rast("OutData/Dp_int")

tiff(filename = 'Maps/Dp_int.tif', units = "in",
     width = ncol(r)/300, 
     height = (ncol(r)/300)/ar, 
     type = "cairo", res = 300, 
     compression = "zip")
par(mai = c(0,0,0,0), omi = c(0,0,0,0))
plot(r, type = "continuous", axes = FALSE, 
     col = turbo(256, direction = -1),
     maxcell = ncell(r), 
     mar = c(2,0,2,4), main = "Potential Diversity Index (Int)")
plot(countries, lwd = 0.4, add = T)
dev.off()




# SDM Models ####################################

## Data prep ---------
# Upload allocated SDM data
d <- fread("OutData/allocated_sdm.csv") 

acols <- names(d)[names(d) %like% '^a\\.']
scols <- gsub("^a", "^s", acols)
crops <- gsub("^a\\.", "", acols)


# total allocated area
d[, tcl:= Reduce(`+`, .SD), .SDcols = acols]

# allocated area fraction
d[, (acols):= .SD/tcl, .SDcols = acols]

# set 0 crop area as NA for plotting
for(j in acols) set(d, which(d[[j]] == 0), j, NA)

# set 0 crop suit as NA for plotting
for(j in scols) set(d, which(d[[j]] == 0), j, NA)

# empty raster
r <- rast(totcl)
v <- rep(NA_real_, ncell(r))

## crop area ------------
for(i in 1:length(acols)){
  v[d$cell] <- d[[acols[i]]]
  values(r) <- v
  
  tiff(filename = file.path(dir_area, paste0(crops[i], "_sdm.tif")),  
       units = "in",
       width = ncol(r)/300, 
       height = (ncol(r)/300)/ar, 
       type = "cairo", res = 300, 
       compression = "zip")
  par(mai = c(0,0,0,0), omi = c(0,0,0,0))
  plot(r, type = "continuous", axes = FALSE,
       col = mako(256, direction = -1),
       maxcell = ncell(r), mar = c(2,0,2,6), 
       main = paste(crops[i], "allocated area (sdm)"))
  plot(countries, lwd = 0.4, add = T)
  dev.off()
}

## suitability -----------
for(i in 1:length(scols)){
  v[d$cell] <- d[[scols[i]]]
  values(r) <- v
  
  tiff(filename = file.path(dir_suit, paste0(crops[i], "_sdm.tif")),  
       units = "in",
       width = ncol(r)/300, 
       height = (ncol(r)/300)/ar, 
       type = "cairo", res = 300, 
       compression = "zip")
  par(mai = c(0,0,0,0), omi = c(0,0,0,0))
  plot(r, type = "continuous", axes = FALSE,
       col = mako(256, direction = -1),
       maxcell = ncell(r), mar = c(2,0,2,6), 
       main = paste(crops[i], "suitability (sdm)"))
  plot(countries, lwd = 0.4, add = T)
  dev.off()
}

## Diversity -------------
r <- rast("OutData/Dp_sdm")

tiff(filename = 'Maps/Dp_sdm.tif', units = "in",
     width = ncol(r)/300, 
     height = (ncol(r)/300)/ar, 
     type = "cairo", res = 300, 
     compression = "zip")
par(mai = c(0,0,0,0), omi = c(0,0,0,0))
plot(r, type = "continuous", axes = FALSE, 
     col = turbo(256, direction = -1),
     maxcell = ncell(r), 
     mar = c(2,0,2,4), main = "Potential Diversity Index (SDM)")
plot(countries, lwd = 0.4, add = T)
dev.off()



# Ecocrop Model ####################################

## Data prep ---------
# Upload allocated SDM data
d <- fread("OutData/allocated_eco.csv") 

acols <- names(d)[names(d) %like% '^a\\.']
scols <- gsub("^a", "^s", acols)
crops <- gsub("^a\\.", "", acols)


# total allocated area
d[, tcl:= Reduce(`+`, .SD), .SDcols = acols]

# allocated area fraction
d[, (acols):= .SD/tcl, .SDcols = acols]

# set 0 crop area as NA for plotting
for(j in acols) set(d, which(d[[j]] == 0), j, NA)

# set 0 crop suit as NA for plotting
for(j in scols) set(d, which(d[[j]] == 0), j, NA)

# empty raster
r <- rast(totcl)
v <- rep(NA_real_, ncell(r))

## crop area ------------
for(i in 1:length(acols)){
  v[d$cell] <- d[[acols[i]]]
  values(r) <- v
  
  tiff(filename = file.path(dir_area, paste0(crops[i], "_eco.tif")),  
       units = "in",
       width = ncol(r)/300, 
       height = (ncol(r)/300)/ar, 
       type = "cairo", res = 300, 
       compression = "zip")
  par(mai = c(0,0,0,0), omi = c(0,0,0,0))
  plot(r, type = "continuous", axes = FALSE,
       col = mako(256, direction = -1),
       maxcell = ncell(r), mar = c(2,0,2,6), 
       main = paste(crops[i], "allocated area (eco)"))
  plot(countries, lwd = 0.4, add = T)
  dev.off()
}

## suitability -----------
for(i in 1:length(scols)){
  v[d$cell] <- d[[scols[i]]]
  values(r) <- v
  
  tiff(filename = file.path(dir_suit, paste0(crops[i], "_eco.tif")),  
       units = "in",
       width = ncol(r)/300, 
       height = (ncol(r)/300)/ar, 
       type = "cairo", res = 300, 
       compression = "zip")
  par(mai = c(0,0,0,0), omi = c(0,0,0,0))
  plot(r, type = "continuous", axes = FALSE,
       col = mako(256, direction = -1),
       maxcell = ncell(r), mar = c(2,0,2,6), 
       main = paste(crops[i], "suitability (eco)"))
  plot(countries, lwd = 0.4, add = T)
  dev.off()
}

## Diversity -------------
r <- rast("OutData/Dp_eco")

tiff(filename = 'Maps/Dp_eco.tif', units = "in",
     width = ncol(r)/300, 
     height = (ncol(r)/300)/ar, 
     type = "cairo", res = 300, 
     compression = "zip")
par(mai = c(0,0,0,0), omi = c(0,0,0,0))
plot(r, type = "continuous", axes = FALSE, 
     col = turbo(256, direction = -1),
     maxcell = ncell(r), 
     mar = c(2,0,2,4), main = "Potential Diversity Index (Ecocrop)")
plot(countries, lwd = 0.4, add = T)
dev.off()


