library(magrittr)
library(data.table)
library(terra)
library(maxnet)
library(ranger)
library(xgboost)

if(system('hostname', TRUE) == "ESP-RH-9891"){
  setwd("D:/globcropdiv/")
} else if(system('hostname', TRUE) == "LAPTOP-ST129J47"){ 
  setwd("G:/My Drive/globcropdiv/")
} # else if { ... 

# functions ----
source("Functions/maxent_output_transformation.R")
source("Functions/sample_presences.R")

# DATA PREP ###########################################################
## Total Cropland -----
d <-"InData/TotalCropland.tif" %>% 
  rast() %>% 
  values() %>% 
  {data.table(cell = which(!is.na(.)), tot = .[!is.na(.)])}

## Crop abundances -----
r_abun <- "InData/CropAbundance/*.tif" %>% Sys.glob() %>% rast() 
crops <- names(r_abun)
d[, (crops):= extract(r_abun, cell)]

# change abundances to fraction
d[, (crops):= lapply(.SD, function(x) x/tot), .SDcols = (crops)]

# crops by source
sc <- fread("AuxData/CropCategories.csv")
monf <- crops[crops %in% sc[source == 'Monf', crop]]
spam <- crops[crops %in% sc[source == 'SPAM', crop]]


## DQI  ------------
# it only affects Monfreda Crops. 
rdqi <- paste0("InData/CropAbundance/DataQualityIndex/", monf, "_DQI.tif") %>% 
  rast()
monf.dqi <- paste0(monf, "_dqi")
d[, (monf.dqi):= extract(rdqi, cell)]

envCV_to_dqi <- function(x) ifelse(is.na(x), 0, pmax(1 - x/4,0))
for(j in monf.dqi) set(d, j = j, value = envCV_to_dqi(d[[j]]))


## Predictors --------
r_preds <- fread("AuxData/SelectedPreds.csv")$file_name %>% 
  file.path("InData", .) %>% 
  Sys.glob() %>% 
  rast()
preds <- fread("AuxData/SelectedPreds.csv")$short_name
d[, (preds):= extract(r_preds, cell)]

# output dir 
outdir <- "OutData/SDM"


# RUN MODELS  ############################################################

## Maxent  --------
me_tp <- readRDS("OutData/SDM/CV/Tuning/ME_selected_pars")
ss <- me_tp$ss
regmult <- me_tp$regmult
npr <- me_tp$npr

tr <- which(d$cell%%(1/ss) == 0)

set.seed(1)

for(i in 1:length(crops)){
  # run model
  pm <- if(crops[i] %in% monf) 3 else 2 
  pr <- fpr(d[tr], crops[i], pm, npr)
  m <- maxnet(c(rep(1, npr), rep(0, d[tr][,.N])),
              d[tr][c(pr, 1:length(tr)), ..preds],
              regmult = regmult)  
  logp <- predict(m, d[, ..preds], type = "logistic")
  
  # write predictions on disc
  r <- rast(r_abun[[1]])
  v <- rep(NA_real_, ncell(r))
  v[d$cell] <- logp
  values(r) <- v
  fn <- file.path(outdir, paste0(crops[i], "_ME.tif"))
  wopts <- list(names = crops[i], filetype = "GTiff",
                gdal=c("COMPRESS=Deflate","PREDICTOR=1", "ZLEVEL=6"))
  writeRaster(r, filename = fn, overwrite = T, wopt = wopts)
  
  rm(pr, m, logp, r, v, fn, wopts); gc()
}


## Random Forest --------
rf_tp <- readRDS("OutData/SDM/CV/Tuning/RF_selected_pars")
mtry <- rf_tp$mtry
min.node.size <- rf_tp$min.node.size

set.seed(0)
for(i in 1:length(crops)){  
  rf <- ranger(x = d[, ..preds], 
               y = d[[crops[i]]],
               num.trees = 1000,
               mtry = mtry,
               min.node.size = min.node.size, 
               num.threads = 4)
  pv <- predict(rf, d[, ..preds])$predictions
  
  # write predictions on disc
  r <- rast(r_abun[[1]])
  v <- rep(NA_real_, ncell(r))
  v[d$cell] <- pv
  values(r) <- v
  fn <- file.path(outdir, paste0(crops[i], "_RF.tif"))
  wopts <- list(names = crops[i], filetype = "GTiff",
                gdal=c("COMPRESS=Deflate","PREDICTOR=1", "ZLEVEL=6"))
  writeRaster(r, filename = fn, overwrite = T, wopt = wopts)
  rm(rf, pv, r, v, fn, wopts); gc()
}


## Boosted Regression Trees -----
brt_tp <- readRDS("OutData/SDM/CV/Tuning/BRT_selected_pars")
brt_tp$nthread <- 4
nrounds <- brt_tp$nrounds
brt_tp$nrounds <- NULL

dm <- as.matrix(d[, ..preds])

set.seed(1)
for(i in 1:length(crops)){  
  brt <- xgboost(params = brt_tp,
                 data = dm,
                 label = d[[crops[i]]],
                 nrounds = nrounds, 
                 verbose = 0)
  pv <- predict(brt, dm)
  # limit output
  pv[pv < 0] <- 0
  pv[pv > 1] <- 1
   
  # write predictions on disc
  r <- rast(r_abun[[1]])
  v <- rep(NA_real_, ncell(r))
  v[d$cell] <- pv
  values(r) <- v
  fn <- file.path(outdir, paste0(crops[i], "_BRT.tif"))
  wopts <- list(names = crops[i], filetype = "GTiff",
                gdal=c("COMPRESS=Deflate","PREDICTOR=1", "ZLEVEL=6"))
  writeRaster(r, filename = fn, overwrite = T, wopt = wopts)
  rm(brt, pv, r, v, fn, wopts); gc()
}
