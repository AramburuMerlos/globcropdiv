library(magrittr)
library(data.table)
library(terra)
library(Recocrop)

if(system('hostname', TRUE) == "ESP-RH-9891"){
  setwd("D:/globcropdiv/")
} else if(system('hostname', TRUE) == "LAPTOP-ST129J47"){ 
  setwd("G:/My Drive/globcropdiv/")
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
dpars[, SID:= paste0("SID_", SID)]
crops_sp <- unique(dpars$SID)

rpar <- c("RMIN", "ROPMN", "ROPMX", "RMAX")
phpar <- c("PHMIN", "PHOPMN", "PHOPMX", "PHMAX")
tpar <- c("TMIN","TOPMN","TOPMX","TMAX")
mrpar <- paste0(rpar, "_M")


# RUN model -----
for(i in 1:length(crops_sp)){
  # rows with parameters for all species in crop[i]
  ssp <- which(dpars$SID == crops_sp[i])
  score_ssp <- paste0("score_", ssp) 
  
  # run for each species within crop category
  for(j in 1:length(ssp)){
    # crop subspecies parameters
    crop_pars <- cbind(
      duration = c(dpars[ssp[j], duration], NA, NA, NA),
      ktmp = dpars[ssp[j], KTMP] + c(-1, +1, NA, NA),
      tavg = unlist(dpars[ssp[j], ..tpar], use.names = FALSE),
      prec = unlist(dpars[ssp[j], ..mrpar], use.names = FALSE),
      ph = unlist(dpars[ssp[j], ..phpar], use.names = FALSE),
      anpr = unlist(dpars[ssp[j], ..rpar], use.names = FALSE)
    )
    
    foo <- list(name = dpars[ssp[j], NAME], 
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
    
    # subspecies prediction
    d[, (score_ssp[j]):= rfp * (1 - irri/100) + irp * irri/100]
    d[, c("rfp", "irp"):= NULL]
  }
  # get maximum score for all subspecies within crop species i
  d[, (crops_sp[i]):= do.call(pmax, .SD), .SDcols = score_ssp]
  d[, (score_ssp):= NULL]
  print(paste(i, "/", length(crops_sp)))
}


d[, rs:= Reduce(`+`, .SD), .SDcols = (crops_sp)]
d[, (crops_sp):= .SD/rs, .SDcols = (crops_sp)]

# diversity function 
fd <- function(x, na.rm = na.rm){
  if(na.rm == TRUE){
    x <- x[x > 0 & !is.na(x)]  
  } else {
    x <- x[x > 0]
  }
  return(exp(-sum(x * log(x))))
}

d[, Dp:= apply(.SD, 1, fd, na.rm = TRUE), .SDcols = crops_sp]
d[, summary(Dp)]

r <- rast(totcl)
v <- rep(NA_real_, ncell(totcl))
v[d$cell] <- d$Dp
values(r) <- v
plot(r)

writeRaster(r, filename = "OutData/pot_D_eco_sp.tif", overwrite = T, 
            wopt = list(names = "pot_D_eco_sp", filetype = "GTiff",
                        gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6"))
)

