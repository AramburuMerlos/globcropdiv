library(magrittr)
library(data.table)
library(terra)
library(maxnet)
library(ranger)
library(xgboost)

if(system('hostname', TRUE) %in% c("ESP-RH-9891", "LAPTOP-ST129J47")){
  setwd("D:/globcropdiv/")
} # else if { ... 


# functions ----
source("Functions/maxent_output_transformation.R")
source("Functions/sample_presences.R")

frmse <- function(obs, pred) sqrt(mean((obs-pred)^2, na.rm = T))

# DATA PREP ###########################################################
## Predictors --------------
preds <- fread("AuxData/SelectedPreds.csv")$file_name %>% 
  file.path("InData", .) %>% 
  Sys.glob() %>% 
  rast()
names(preds) <- fread("AuxData/SelectedPreds.csv")$short_name

## Total Cropland (per cell) --------
# predictions will be made for cropland only
totcl <- rast("InData/TotalCropland.tif")

## Set train/test regions -----------
# define regions extent
cext <- as.vector(ext(totcl))                     # Roughly:
lext <- list(ext(cext[1],-30 , 10, cext[4]),      # North America
             ext(cext[1],-30 , cext[3], 10),      # South America
             ext(-30, 60, cext[3], 20),           # Sub-Sahara Africa  
             ext(-30, 60, 20, cext[4]),           # N.Africa + Europe + E.Asia
             ext(60, cext[2], 5, cext[4]),        # M&E (N) Asia
             ext(60, cext[2], cext[3], 5))        # SE.Asia + Australia/Oceania

nk <- length(lext)
# get cell number and total cropland values for each region, create DT 
ld <- vector(mode = "list", length = nk)
for (i in seq(nk)){
  ld[[i]] <- data.table(reg = i, cell = cells(totcl, lext[[i]]))
  ld[[i]][, tot:= extract(totcl, cell)]
} 
d <- rbindlist(ld)
rm(ld)

# remove NA (non cropland) cells, order by cell number
d <- d[!is.na(tot),]
setorder(d, cell)


## Select Crop Species for testing------
# the model will be tested on evenly widespread crops 

# Crop abundance data
r_abun <- "InData/CropAbundance/*.tif" %>% Sys.glob() %>% rast() 
crops <- names(r_abun)
d[, (crops):= extract(r_abun, cell)]

# get total abundance for each crop at each region
regtots <- d[, lapply(.SD, sum), .SDcols = (crops), by = reg]
regtots <- melt(regtots, id.vars = "reg",  measure.vars = crops, 
                variable.name = "crop", value.name = "area") 

# crop area fraction of total cropland of each region
regtots[, frac:= area/sum(area), by = reg]

# coefficient of variation and minimum crop fraction across regions
regvar <- regtots[, .(CV = sd(frac)/mean(frac),
                      min_frac = min(frac)), by = crop]
                     
# select testing crops
crops <- regvar[CV < .8 & min_frac > 0.001, crop] 
crops <- as.character(crops)

# remove non-selected crops from d and rabun
ctorm <- names(r_abun)[!names(r_abun) %in% crops]
d[, (ctorm):= NULL]
r_abun <- r_abun[crops]

# crops by source
sc <- fread("AuxData/CropCategories.csv")
monf <- crops[crops %in% sc[source == 'Monf', crop]]
spam <- crops[crops %in% sc[source == 'SPAM', crop]]

# change abundances to fraction
d[, (crops):= lapply(.SD, function(x) x/tot), .SDcols = (crops)]

## Extract predictor values --------
preds_names <- names(preds)
d[, (preds_names):= extract(preds, cell)]


## get DQI  ------------
# it only affects Monfreda Crops. 
rdqi <- paste0("InData/CropAbundance/DataQualityIndex/", monf, "_DQI.tif") %>% 
 rast()
monf.dqi <- paste0(monf, "_dqi")

## Extract DQI values --------
d[, (monf.dqi):= extract(rdqi, cell)]

# DQI is a measure of environmental variability. Must be changed to an index
# NA means either that there was no Monfreda data or that is outside used country borders
envCV_to_dqi <- function(x) ifelse(is.na(x), 0, pmax(1 - x/4,0))
for(j in monf.dqi) set(d, j = j, value = envCV_to_dqi(d[[j]]))

# output dir for CV results
tundir <- "OutData/SDM/CV/Tuning"
dir.create(tundir, F, T)


# TUNE MODELS #################################################################
## 0. NULL MODEL ------
# average crop fraction 
null_rmse <- function(d, crops){
  rmse <- data.table(crops = crops)
  rmse[, null:= mapply(frmse, 
                       d[,.SD, .SDcols = (crops)], 
                       d[,lapply(.SD,mean), .SDcols = (crops)])]
}

## 1. MAXENT ------------

### i. Presence method ----
rmse <- null_rmse(d, crops)
cors <- data.table(crops = crops)

# presence methods
pms <- 1:5

# 2% of the data will be used for model training. 
# same cells for all crops and models. Sampled regularly
ss = 0.02
tr <- which(d$cell%%(1/ss) == 0)

# default parameters:
nbg <- 8000
npr <- 2000
regmult <- 2

set.seed(0)
for(i in 1:length(crops)){
  for(pm in pms){
    if(pm %in% c(3,5) & !(crops[i] %in% monf)){
      next
    } 
    for(k in 1:nk){
      pr <- fpr(d[tr], crops[i], pm, npr, colk = 'reg', k = k)
      bg <- sample(d[tr][, .I[reg != k]], nbg)
      m <- maxnet(c(rep(1, npr), rep(0, nbg)),
                  d[tr][c(pr, bg), ..preds_names],
                  regmult = regmult)  
      lp <- predict(m, d[reg == k, ..preds_names], type = "link")
      mean_abun <- mean(d[[crops[i]]][d$reg == k])
      clog <- clogtr(mean_abun, lp)
      d[reg == k, paste0(crops[i],".ME"):= clog]
      rm(pr, m, lp, clog, mean_abun); gc()
      print(paste0("ME: ", i, ".", pm, ".",  k))
    }
    set(rmse, i, paste0("ME.pm", pm),
        frmse(d[[crops[i]]], d[[paste0(crops[i],".ME")]]))
    set(cors, i, paste0("ME.pm", pm),
        cor(d[[crops[i]]], d[[paste0(crops[i],".ME")]]))
  }
  d[, paste0(crops[i], ".ME"):= NULL] #clean up predictions
}

fwrite(rmse, file.path(tundir, "ME_pm_rmse.csv"))
fwrite(cors, file.path(tundir, "ME_pm_cors.csv"))

colMeans(rmse[,-1]) %>% which.min() %>% names()
## [1] "ME.pm2"
colMeans(cors[,-1]) %>% which.max() %>% names()
## [1] "ME.pm2"

# check if DQI improves correlation
cors[!is.na(ME.pm3), ME.pm3 > ME.pm2]
## [1] TRUE


### ii. Tuning parameters -----
# Sampling regularly and then sampling presences from that sample might be better than sampling presences from all data. Maybe because this generates presences that better explore all environments. 

rmse <- null_rmse(d, crops)
cors <- data.table(crops = crops)


# tuning parameters
nbg <- 8000
tp <- expand.grid(regmult = 2^(0:3),                # shrinkage
                  ss = c(0.02, 0.1, 0.5, 1),        # regular sampling
                  npr = c(1000, 2000, 4000))        # number of presences

set.seed(0)
for(i in 1:length(crops)){
  pm <- if(crops[i] %in% monf) 3 else 2 
  for(j in 1:nrow(tp)){
    
    # is this one already done?
    if(!is.na(rmse[[paste0("ME.tp", j)]][i])) next

    for(k in 1:nk){
      tr <- which(d$cell %% (1/tp$ss[j]) == 0)
      pr <- fpr(d[tr], crops[i], pm, tp$npr[j], colk = 'reg', k = k)
      bg <- d[tr][, sample(.I[reg != k], nbg, replace = sum(reg != k) < nbg)]
      m <- maxnet(c(rep(1, tp$npr[j]), rep(0, nbg)),
                  d[tr][c(pr, bg), ..preds_names],
                  regmult = tp$regmult[j])  
      lp <- predict(m, d[reg == k, ..preds_names], type = "link")
      mean_abun <- mean(d[[crops[i]]][d$reg == k])
      clog <- clogtr(mean_abun, lp)
      d[reg == k, paste0(crops[i],".ME"):= clog]
      rm(pr, m, lp, clog, mean_abun); gc()
      print(paste0("ME: ", i, ".", j, ".",  k))
    }
    set(rmse, i, paste0("ME.tp", j),
        frmse(d[[crops[i]]], d[[paste0(crops[i],".ME")]]))
    set(cors, i, paste0("ME.tp", j),
        cor(d[[crops[i]]], d[[paste0(crops[i],".ME")]]))
  }
  d[, paste0(crops[i], ".ME"):= NULL] #clean up predictions
}

fwrite(rmse, file.path(tundir, "ME_tp_rmse.csv"))
fwrite(cors, file.path(tundir, "ME_tp_cors.csv"))
fwrite(tp, file.path(tundir, "ME_tunpar.csv"))



## 2. RF ---------------
# 2% of the data will be used for model training. 
# same cells for all crops and models. Sampled regularly
ss = 0.02
tr <- which(d$cell%%(1/ss) == 0)

rmse <- null_rmse(d, crops)
cors <- data.table(crops = crops)

tp <- expand.grid(
  mtry = c(1:4,seq(6,nlyr(preds),4)),
  min.node.size = seq(5,50,5)
)

set.seed(0)
for(i in 1:length(crops)){  
  for(j in 1:nrow(tp)){
    for(k in 1:nk){
      rf <- ranger(x = d[tr][reg != k, ..preds_names], 
                   y = d[tr][reg != k][[crops[i]]],
                   num.trees = 500,
                   mtry = tp$mtry[j],
                   min.node.size = tp$min.node.size[j])
      pv <- predict(rf, d[reg == k, ..preds_names])$predictions
      d[reg == k, paste0(crops[i],".RF"):= pv]
      print(paste0("RF: ", i, ".", j, ".",  k))
    }
    set(rmse, i, paste0("RF.tp.", j),
        frmse(d[[crops[i]]], d[[paste0(crops[i],".RF")]]))
    set(cors, i, paste0("RF.tp.", j),
        cor(d[[crops[i]]], d[[paste0(crops[i],".RF")]]))
  }
  d[, paste0(crops[i], ".RF"):= NULL] #clean up predictions
}

fwrite(rmse, file.path(tundir, "RF_rmse.csv"))
fwrite(cors, file.path(tundir, "RF_cors.csv"))
fwrite(tp, file.path(tundir, "RF_tunpar.csv"))


## 3. BRT ---------------
# 2% of the data will be used for model training. 
# same cells for all crops and models. Sampled regularly
ss = 0.02
tr <- which(d$cell%%(1/ss) == 0)

rmse <- null_rmse(d, crops)
cors <- data.table(crops = crops)

tp <- expand.grid(
  eta = c(.01, 0.02, .1),
  max_depth = c(1:3, 5),
  min_child_weight = c(10, 50, 100),
  subsample = c(.5, 1), 
  colsample_bytree = c(.2, .4, .8),
  nrounds = seq(400,1200,400)              
)

set.seed(0)
for(i in 1:length(crops)){  
  for(j in 1:nrow(tp)){
    for(k in 1:nk){
      params <- list(eta = tp$eta[j],
                     max_depth = tp$max_depth[j],
                     min_child_weight = tp$min_child_weight[j],
                     subsample = tp$subsample[j],
                     colsample_bytree = tp$colsample_bytree[j], 
                     nthread = 4
      )
      brt <- xgboost(params = params,
                     data = as.matrix(d[tr][reg != k, ..preds_names]), 
                     label = d[tr][reg != k][[crops[i]]],
                     nrounds = tp$nrounds[j], 
                     verbose = 0)
      pv <- predict(brt, as.matrix(d[reg == k, ..preds_names]))
      d[reg == k, paste0(crops[i],".BRT"):= pv]
      print(paste0("BRT: ", i, ".", j, ".",  k))
    }
    set(rmse, i, paste0("BRT.tp.", j),
        frmse(d[[crops[i]]], d[[paste0(crops[i],".BRT")]]))
    set(cors, i, paste0("BRT.tp.", j),
        cor(d[[crops[i]]], d[[paste0(crops[i],".BRT")]]))
  }
  d[, paste0(crops[i], ".BRT"):= NULL] #clean up predictions
}

fwrite(rmse, file.path(tundir, "BRT_rmse.csv"))
fwrite(cors, file.path(tundir, "BRT_cors.csv"))
fwrite(tp, file.path(tundir, "BRT_tunpar.csv"))


# COMPARE MODELS  ############################################################

## upload results ------
tundir <- "OutData/SDM/CV/Tuning"

rmse_ME <- fread(file.path(tundir, "ME_tp_rmse.csv"))
cors_ME <- fread(file.path(tundir, "ME_tp_cors.csv"))

rmse_RF <- fread(file.path(tundir, "RF_rmse.csv"))
cors_RF <- fread(file.path(tundir, "RF_cors.csv"))

rmse_BRT <- fread(file.path(tundir, "BRT_rmse.csv"))
cors_BRT <- fread(file.path(tundir, "BRT_cors.csv"))

tp_ME <- fread(file.path(tundir, "ME_tunpar.csv"))
tp_RF <- fread(file.path(tundir, "RF_tunpar.csv"))
tp_BRT <- fread(file.path(tundir, "BRT_tunpar.csv"))

## winners ----

# which method had the best accuracy?
# reduction of null rmse (percentage)
f <- function(x, null) (null - x)/ null * 100

# RMSE reduction from NULL model (%)
list(rmse_ME[,-(1:2)], rmse_RF[,-(1:2)], rmse_BRT[,-(1:2)]) %>% 
  sapply(colMeans) %>% 
  sapply(min) %>% 
  f(mean(rmse_ME$null))

# correlation 
list(cors_ME[,-(1:2)], cors_RF[,-(1:2)], cors_BRT[,-(1:2)]) %>% 
  sapply(colMeans) %>% 
  sapply(max)

### Maxent has the highest correlation and RMSE reduction

# parameter combinations with lowest RMSE and greatest correlation
ME_min_rmse <- rmse_ME[,-1] %>% colMeans() %>% which.min() %>% names()
RF_min_rmse <- rmse_RF[,-1] %>% colMeans() %>% which.min() %>% names()
BRT_min_rmse <- rmse_BRT[,-1] %>% colMeans() %>% which.min() %>% names()

ME_max_cors <- cors_ME[,-1] %>% colMeans() %>% which.max() %>% names()
RF_max_cors <- cors_RF[,-1] %>% colMeans() %>% which.max() %>% names()
BRT_max_cors <- cors_BRT[,-1] %>% colMeans() %>% which.max() %>% names()

tp_ME[as.numeric(gsub("ME.tp", "", ME_min_rmse)),]
tp_ME[as.numeric(gsub("ME.tp", "", ME_max_cors)),]
# same ss and npr; different regmult

tp_RF[as.numeric(gsub("RF.tp.", "", RF_min_rmse)),]
tp_RF[as.numeric(gsub("RF.tp.", "", RF_max_cors)),]
# same mtry, different min.node.size

tp_BRT[as.numeric(gsub("BRT.tp.", "", BRT_min_rmse)),]
tp_BRT[as.numeric(gsub("BRT.tp.", "", BRT_max_cors)),]
# same eta, max_depth, colsample_bytree and nrounds; different min_child_weight and subsample

# Chose one or average differences? Is there a clear winner?
# how do they compare at the accuracy measurement that they did worse?


### Maxent  ---- 
colMeans(cors_ME[, ..ME_min_rmse])/ max(colMeans(cors_ME[, -1]))

f(colMeans(rmse_ME[, ..ME_max_cors]), mean(rmse_ME$null)) / f(min(colMeans(rmse_ME[, -1])), mean(rmse_ME$null)) 

# almost the same difference, use average regmult
rbind(
  tp_ME[as.numeric(gsub("ME.tp", "", ME_min_rmse)),],
  tp_ME[as.numeric(gsub("ME.tp", "", ME_max_cors)),]
) %>% 
  as.list() %>% 
  lapply(mean) %>% 
  saveRDS(file.path(tundir,"ME_selected_pars"))


### Random Forest ----
colMeans(cors_RF[, ..RF_min_rmse])/ max(colMeans(cors_RF[, -1]))

f(colMeans(rmse_RF[, ..RF_max_cors]), mean(rmse_RF$null)) / f(min(colMeans(rmse_RF[, -1])), mean(rmse_RF$null)) 

# the parameters that gave the lowest RMSE are the winners
tp_RF[as.numeric(gsub("RF.tp.", "", RF_min_rmse)),] %>% 
  as.list() %>% 
  saveRDS(file.path(tundir,"RF_selected_pars"))


### Boosted Regression Trees ----
colMeans(cors_BRT[, ..BRT_min_rmse])/ max(colMeans(cors_BRT[, -1]))

f(colMeans(rmse_BRT[, ..BRT_max_cors]), mean(rmse_BRT$null)) / f(min(colMeans(rmse_BRT[, -1])), mean(rmse_BRT$null)) 

# the parameters that gave the lowest RMSE are the winners
tp_BRT[as.numeric(gsub("BRT.tp.", "", BRT_min_rmse)),]%>% 
  as.list() %>% 
  saveRDS(file.path(tundir,"BRT_selected_pars"))

### Averaging weights -------
# Based on rmse reduction from the null
w <- c(colMeans(rmse_ME[, ..ME_min_rmse]),
       colMeans(rmse_RF[, ..RF_min_rmse]),
       colMeans(rmse_BRT[, ..BRT_min_rmse])) %>% 
  f(mean(rmse_ME$null)) %>% 
  `/`(sum(.)) %>% 
  unname()

names(w) <- c("ME", "RF", "BRT")
saveRDS(w, file.path(tundir, "Model_averaging_weights.RDS"))
