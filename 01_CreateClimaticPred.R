library(raster)
library(terra)
library(data.table)

nwd <- nchar(getwd())
if(substr(getwd(),  nwd - 10, nwd) != "globcropdiv") stop("See 0000_wd.R")

if(dir.exists("D:/WorldClim")){
  wcpath <- "D:/WorldClim/2.1/wc5min"
} else {
  wcpath <- "InData/WorldClim/2.1/wc5min"
}

outpath <- file.path(wcpath, "extra")
dir.create(outpath)

# Growing Degree Days ---------
tavg <- rast(Sys.glob(file.path(wcpath, "tavg/*.tif")))
days_per_month <- c(31,28.25,31,30,31,30,31,31,30,31,30,31) 
fgdd <- function(x){
  x[x < 0] <- NA
  return(sum(x * days_per_month, na.rm = T))
}

app(tavg, fgdd, 
    filename = file.path(outpath, "GDD.tif"), 
    overwrite = T, 
    wopt = list(names = "GDD", filetype = "GTiff", 
                gdal = c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))

# Radiation on top of the atmosphere (monthly average) -------------
## latitude raster
lat <- yFromCell(raster(tavg), 1:ncell(tavg))
rlat <- tavg[[1]]
rlat <- setValues(rlat, lat)
rlat <- mask(rlat, tavg[[1]])
rm(lat)

# function to estimate average top of atmosphere solar radiation in MJ/day
source("Functions/solarrad.R")

for (i in 1:12){
  app(rlat, rta, month = i, nodes = 8,      
      filename = file.path(outpath, paste0("RadTopAtm_m", formatC(i, width = 2, flag = '0'), ".tif")), 
      overwrite = T, 
      wopt = list(names = paste0("rta", i), filetype = "GTiff", progress = 1, 
                  gdal = c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6"))) 
}

# Potential Evapotranspiration (monthly average) ------------
rta <- rast(Sys.glob(file.path(outpath, "RadTopAtm_m*.tif"))) # in MJ/day
rta_mm <- rta/2.45 # mm/day (see http://www.fao.org/3/X0490E/x0490e07.htm)
tavg <- rast(Sys.glob(file.path(wcpath, "tavg/*.tif")))
tmin <- rast(Sys.glob(file.path(wcpath, "tmin/*.tif")))
tmax <- rast(Sys.glob(file.path(wcpath, "tmax/*.tif")))
tran <- tmax - tmin

fpet <- function(rta_mm, tavg, tran) 0.0023 * rta_mm * (tavg + 17.8) * sqrt(tran) # [mm/day]

for(i in 1:12){
  lapp(c(rta_mm[[i]], tavg[[i]], tran[[i]]), fun = fpet, 
       filename = file.path(outpath, paste0("PET_m", formatC(i, width = 2, flag = '0'), ".tif")), 
       overwrite = T, 
       wopt = list(names = paste0("rta", i), filetype = "GTiff",
                   gdal = c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6"))) 
}

# Aridity Index (annual) ------------
pet <- rast(Sys.glob(file.path(outpath, "PET_m*.tif"))) # daily average
days_per_month <- c(31,28.25,31,30,31,30,31,31,30,31,30,31) 
pet <- pet * days_per_month # mm per month
prec <- rast(Sys.glob(file.path(wcpath, "prec/*.tif"))) # mm per month

pet_annual <- sum(pet)
prec_annual <- sum(prec)

lapp(c(prec_annual, pet_annual), fun = function(x,y)pmin(pmax(x/y,0),10), 
     filename = file.path(outpath, "AI.tif"), 
     overwrite = T, 
     wopt = list(names = "Aridity_Index", filetype = "GTiff",
                 gdal = c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6"))) 

ai <- rast(file.path(outpath, "AI.tif"))
plot(ai)
# Pet seasonality ------------
fps <- function(x) sd(x, na.rm = T)
app(pet, fps, 
    filename = file.path(outpath, "PET_seasonality.tif"), 
    overwrite = T, 
    wopt = list(names = "PET_seasonality", filetype = "GTiff",
                gdal = c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))
ps <- rast(file.path(outpath, "PET_seasonality.tif"))

