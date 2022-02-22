library(magrittr)
library(data.table)
library(terra)
library(Recocrop)

if(system('hostname', TRUE) %in% c("ESP-RH-9891", "LAPTOP-ST129J47")){
  setwd("D:/globcropdiv/")
} # else if { ... 

wcdir <- "InData/WorldClim/2.1/wc5min"
sdir <- "InData/SoilGrids"
outdir <- "OutData/Ecocrop"
dir.create(outdir, F, T)

# Load Data -----------
# Prepare data table
totcl <- rast("InData/TotalCropland.tif")
d <- data.table(cell = which(!is.na(values(totcl))))

# Minimum (Killing) Temperature
r_tmin <- "tmin/*.tif" %>% file.path(wcdir, .) %>% Sys.glob() %>% rast()
tmin <- paste0("tmin_", 1:12) 
d[, (tmin):= extract(r_tmin, cell)]

# Monthly Avg temp
r_tavg <- "tavg/*.tif" %>% file.path(wcdir, .) %>% Sys.glob() %>% rast()
tavg <- paste0("tavg_", 1:12) 
d[, (tavg):= extract(r_tavg, cell)]

# Precipitation
r_prec <- "prec/*.tif" %>% file.path(wcdir, .) %>% Sys.glob() %>% rast()
prec <- paste0("prec_", 1:12) 
d[, (prec):= extract(r_prec, cell)]

# Annual precipitation
d[, anpr:= Reduce(`+`, .SD), .SDcols = prec]

# irrigation
r_irri <- rast("InData/AQUASTAT/gmia_v5_aei_pct.asc")
d[, irri:= extract(r_irri, cell)]

# ph
r_ph <- rast(file.path(sdir, "phh2o/phh2o_0-15cm_mean_5min.tif"))
d[, phx10:= extract(r_ph, cell)]
d[, ph:= phx10/10]
d[, phx10:=NULL]

rm(list = ls()[ls() %like% "r_"])

# Parameters --------
dpars <- fread("AuxData/CalibEcoPars.csv")
crops <- unique(dpars$crop)

rpar <- c("RMIN", "ROPMN", "ROPMX", "RMAX")
phpar <- c("PHMIN", "PHOPMN", "PHOPMX", "PHMAX")
tpar <- c("TMIN","TOPMN","TOPMX","TMAX")
mrpar <- paste0(rpar, "_M")


# RUN model -----
for(i in 1:length(crops)){
  # rows with parameters for all species in crop[i]
  spp <- which(dpars$crop == crops[i])
  score_spp <- paste0("score_", spp) 
  
  # run for each species within crop category
  for(j in 1:length(spp)){
    # crop species parameters
    crop_pars <- cbind(
      duration = c(dpars[spp[j], duration], NA, NA, NA),
      ktmp = dpars[spp[j], KTMP] + c(-1, +1, NA, NA),
      tavg = unlist(dpars[spp[j], ..tpar], use.names = FALSE),
      prec = unlist(dpars[spp[j], ..mrpar], use.names = FALSE),
      ph = unlist(dpars[spp[j], ..phpar], use.names = FALSE),
      anpr = unlist(dpars[spp[j], ..rpar], use.names = FALSE)
    )
    
    foo <- list(name = dpars[spp[j], NAME], 
                parameters = crop_pars)
    # rainfed 
    m <- ecocrop(foo)
    control(m, get_max=TRUE)
    
    dynamicPredictors(m) <- cbind(
      ktmp = as.vector(t(d[, ..tmin])),
      tavg = as.vector(t(d[, ..tavg])),
      prec = as.vector(t(d[, ..prec]))
    ) 
    staticPredictors(m) <- cbind(
      ph = d$ph,
      anpr = d$anpr
    )
    d[, rfp:= run(m)]
    
    # irrigated
    m <- ecocrop(foo)
    control(m, get_max=TRUE)
    dynamicPredictors(m) <- cbind(
      ktmp = as.vector(t(d[, ..tmin])),
      tavg = as.vector(t(d[, ..tavg]))
    ) 
    staticPredictors(m) <- cbind(ph = d$ph)
    d[, irp:= run(m)]
    
    # species prediction
    d[, (score_spp[j]):= rfp * (1 - irri/100) + irp * irri/100]
    d[, c("rfp", "irp"):= NULL]
  }
  # get maximum score for all species within crop i
  d[, (crops[i]):= do.call(pmax, .SD), .SDcols = score_spp]
  d[, (score_spp):= NULL]
  
  # write values to spatRast and disk
  r <- rast(totcl)
  v <- rep(NA_real_, ncell(r))
  v[d$cell] <- d[[crops[i]]]
  values(r) <- v
  writeRaster(r, filename = file.path(outdir, paste0(crops[i], ".tif")),
              overwrite = T, 
              wopt = list(names = crops[i], filetype = "GTiff",
                          gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6"))
              )
}




  