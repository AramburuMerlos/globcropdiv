library(raster)
library(terra)
library(data.table)

if(!grepl("globcropdiv$", getwd())){
  if (system('hostname', TRUE) %in% c("ESP-RH-9891", "LAPTOP-ST129J47")) { 
    setwd("G:/My Drive/globcropdiv/")
  } # else if { ... 
}

wcpath <- "InData/WorldClim/2.1/wc5min"

outpath <- file.path(wcpath, "extra")
dir.create(outpath, F, T)

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
      filename = file.path(outpath, 
                           paste0("RadTopAtm_m", 
                                  formatC(i, width = 2, flag = '0'),
                                  ".tif")), 
      overwrite = T, 
      wopt = list(names = paste0("rta", i), 
                  filetype = "GTiff", progress = 1, 
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

# Aridity Index --------------------------------------------------------
pet <- rast(Sys.glob(file.path(outpath, "PET_m*.tif"))) # daily average
days_per_month <- c(31,28.25,31,30,31,30,31,31,30,31,30,31) 
pet <- arith(pet, function(x) x * days_per_month) # mm per month
prec <- rast(Sys.glob(file.path(wcpath, "prec/*.tif"))) # mm per month

# *Annual -------------
pet_annual <- app(pet, sum)
prec_annual <- app(prec, sum)

lapp(c(prec_annual, pet_annual), fun = function(x,y)pmin(pmax(x/y,0),10), 
     filename = file.path(outpath, "AIann.tif"), 
     overwrite = T, 
     wopt = list(names = "AIann", filetype = "GTiff",
                 gdal = c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6"))) 

# *By quarters (3m), starting every month ----
qi <- vector(mode = "list", length = 12)
for(i in 1:12){
  is <- i:(i+2)
  is[is>12] <- is[is>12] - 12 
  qi[[i]] <- is
} 

faiq <- function(x) pmin(pmax(prod(x),0),10)
for (i in 1:12){
  prec_q <- prec[[qi[[i]]]]
  pet_q <- arith(pet[[qi[[i]]]], function(x) 1/x)
  app(tapp(c(prec_q, pet_q), 1:3, fun = faiq), mean, na.rm = T,
      filename = file.path(outpath, 
                           paste0("AIq_q", 
                                  formatC(i, width = 2, flag = '0'),
                                  ".tif")),
      overwrite = T, 
      wopt = list(names = paste0("AIq", i), filetype = "GTiff",
                  gdal = c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))
}


# + AI driest quarter -------
aiq <- rast(Sys.glob(file.path(outpath,"AIq_q*tif")))
app(aiq, min, na.rm = T,
    filename = file.path(outpath, "AIdq.tif"),
    overwrite = T, 
    wopt = list(names = "AIdq", filetype = "GTiff",
                gdal = c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))

# + AI wettest quarter -------
app(aiq, max, na.rm = T,
    filename = file.path(outpath, "AIwq.tif"),
    overwrite = T, 
    wopt = list(names = "AIwq", filetype = "GTiff",
                gdal = c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))

# __(tavg by q) -------
# there might be a better way? but this is very fast
qtavg <- c(
  app(c(tavg[[1]], tavg[[2]], tavg[[3]]), mean, na.rm = T),
  app(c(tavg[[2]], tavg[[3]], tavg[[4]]), mean, na.rm = T),
  app(c(tavg[[3]], tavg[[4]], tavg[[5]]), mean, na.rm = T),
  app(c(tavg[[4]], tavg[[5]], tavg[[6]]), mean, na.rm = T),
  app(c(tavg[[5]], tavg[[6]], tavg[[7]]), mean, na.rm = T),
  app(c(tavg[[6]], tavg[[7]], tavg[[8]]), mean, na.rm = T),
  app(c(tavg[[7]], tavg[[8]], tavg[[9]]), mean, na.rm = T),
  app(c(tavg[[8]], tavg[[9]], tavg[[10]]), mean, na.rm = T),
  app(c(tavg[[9]], tavg[[10]], tavg[[11]]), mean, na.rm = T),
  app(c(tavg[[10]], tavg[[11]], tavg[[12]]), mean, na.rm = T),
  app(c(tavg[[11]], tavg[[12]], tavg[[1]]), mean, na.rm = T),
  app(c(tavg[[12]], tavg[[1]], tavg[[2]]), mean, na.rm = T)
)

# + AI hottest quarter -------
# create matrix to subset AI of hottest (max Tavg) quarter
hqi <- which.max(qtavg)
vhqi <- values(hqi)
submat <- cbind(1:length(vhqi), vhqi)

# get AI by quarter values, create empty raster subset and assign values
vaiq <- values(aiq)
aihq <- rast(hqi)
values(aihq) <- vaiq[submat]

# write on disk
writeRaster(aihq,
            filename = file.path(outpath, "AIhq.tif"),
            overwrite = T, 
            wopt = list(names = "AIhq", filetype = "GTiff",
                        gdal = c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))

# + AI coldest quarter -------
# create matrix to subset AI of coldest (min Tavg) quarter
cqi <- which.min(qtavg)
vcqi <- values(cqi)
submat_c <- cbind(1:length(vcqi), vcqi)

# get AI by quarter values, create empty raster subset and assign values
aicq <- rast(cqi)
values(aicq) <- vaiq[submat_c]

# write on disk
writeRaster(aicq,
            filename = file.path(outpath, "AIcq.tif"),
            overwrite = T, 
            wopt = list(names = "AIcq", filetype = "GTiff",
                        gdal = c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))

# Pet seasonality ------------
fps <- function(x) sd(x, na.rm = T)
app(pet, fps, 
    filename = file.path(outpath, "PET_seasonality.tif"), 
    overwrite = T, 
    wopt = list(names = "PET_seasonality", filetype = "GTiff",
                gdal = c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))
ps <- rast(file.path(outpath, "PET_seasonality.tif"))

