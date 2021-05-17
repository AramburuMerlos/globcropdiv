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




# Diversity -------------

## Actual ----
Da <- rast("OutData/Da.tif")

tiff(filename = 'Maps/Da.tif', units = "in",
     width = ncol(Da)/300, 
     height = (ncol(Da)/300)/ar, 
     type = "cairo", res = 300, 
     compression = "zip")
par(mai = c(0,0,0,0), omi = c(0,0,0,0))
plot(Da, type = "continuous", axes = FALSE, 
     col = turbo(256, direction = -1),
     maxcell = ncell(Da), 
     mar = c(2,0,2,4), main = "Actual Diversity")
plot(countries, lwd = 0.4, add = T)
dev.off()

# Shannon Diversity index for actual diversity (better visualization)
Ha <- log(Da)

tiff(filename = 'Maps/Ha.tif', units = "in",
     width = ncol(Da)/300, 
     height = (ncol(Da)/300)/ar, 
     type = "cairo", res = 300, 
     compression = "zip")
par(mai = c(0,0,0,0), omi = c(0,0,0,0))
plot(Ha, type = "continuous", axes = FALSE, 
     col = turbo(256, direction = -1),
     maxcell = ncell(Da), 
     mar = c(2,0,2,4), main = "Actual Entropy")
plot(countries, lwd = 0.4, add = T)
dev.off()



## Potential -------------
### Integrated ------
Dp_int <- rast("OutData/Dp_int.tif")

tiff(filename = 'Maps/Dp_int.tif', units = "in",
     width = ncol(Dp_int)/300, 
     height = (ncol(Dp_int)/300)/ar, 
     type = "cairo", res = 300, 
     compression = "zip")
par(mai = c(0,0,0,0), omi = c(0,0,0,0))
plot(Dp_int, type = "continuous", axes = FALSE, 
     col = turbo(256, direction = -1),
     maxcell = ncell(Dp_int), 
     mar = c(2,0,2,4), main = "Potential Diversity (Integrated Models)")
plot(countries, lwd = 0.4, add = T)
dev.off()

### SDM  -------------
Dp_sdm <- rast("OutData/Dp_sdm.tif")

tiff(filename = 'Maps/Dp_sdm.tif', units = "in",
     width = ncol(Dp_sdm)/300, 
     height = (ncol(Dp_sdm)/300)/ar, 
     type = "cairo", res = 300, 
     compression = "zip")
par(mai = c(0,0,0,0), omi = c(0,0,0,0))
plot(Dp_sdm, type = "continuous", axes = FALSE, 
     col = turbo(256, direction = -1),
     maxcell = ncell(Dp_sdm), 
     mar = c(2,0,2,4), main = "Potential Diversity (SDM)")
plot(countries, lwd = 0.4, add = T)
dev.off()

### Ecocrop -------------
Dp_eco <- rast("OutData/Dp_eco.tif")

tiff(filename = 'Maps/Dp_eco.tif', units = "in",
     width = ncol(Dp_eco)/300, 
     height = (ncol(Dp_eco)/300)/ar, 
     type = "cairo", res = 300, 
     compression = "zip")
par(mai = c(0,0,0,0), omi = c(0,0,0,0))
plot(Dp_eco, type = "continuous", axes = FALSE, 
     col = turbo(256, direction = -1),
     maxcell = ncell(Dp_eco), 
     mar = c(2,0,2,4), main = "Potential Diversity (Ecocrop)")
plot(countries, lwd = 0.4, add = T)
dev.off()





# Area and Suitability ###################################

## Actual  ----
dir_area <- 'Maps/Area'
dir.create(dir_area, F, T)

dir_suit <- 'Maps/Suit'
dir.create(dir_suit, F, T)

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



## Integrated ----

### Data prep ---------
# Upload allocated SDM data
d <- fread("OutData/allocated_int.csv") 

acols <- names(d)[names(d) %like% '^a\\.']
scols <- gsub("^a", "s", acols)
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

### crop area ------------
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

### suitability -----------
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



## SDM Models ---------

### Data prep ---------
# Upload allocated SDM data
d <- fread("OutData/allocated_sdm.csv") 

acols <- names(d)[names(d) %like% '^a\\.']
scols <- gsub("^a", "s", acols)
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

### crop area ------------
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

### suitability -----------
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




## Ecocrop Model ----------

### Data prep ---------
# Upload allocated SDM data
d <- fread("OutData/allocated_eco.csv") 

acols <- names(d)[names(d) %like% '^a\\.']
scols <- gsub("^a", "s", acols)
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

### crop area ------------
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

### suitability -----------
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


