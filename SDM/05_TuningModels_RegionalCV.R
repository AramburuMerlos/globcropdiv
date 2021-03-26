library(data.table)
library(terra)
library(maxnet)
library(ranger)
library(xgboost)

if(!getwd() %like%  "globcropdiv$") warning("See 0000_wd.R")

# directory of input data
indir <- "D:"

# functions ----
source("Functions/maxent_output_transformation.R")
source("Functions/sample_presences.R")
#source("Functions/row_weight_RF.R")
frmse <- function(obs, pred) sqrt(mean((obs-pred)^2, na.rm = T))
fbrmse <- function(obs, pred, b = 3){
  sqrt(mean(ifelse(obs > pred, 
                   (b * (obs - pred))^2, 
                   (obs - pred)^2), 
            na.rm = T))
}

# DATA PREP ###########################################################
## Predictors --------------
preds.fn <- fread("AuxData/SelectedPreds.csv")$file_name
preds <- rast(file.path(indir, preds.fn))
names(preds) <- fread("AuxData/SelectedPreds.csv")$short_name

## Total Cropland (per cell) --------
# predictions will be made for cropland only
totcl <- rast("InData/TotalCropland.tiff")

# Mask Total Cropland based on predictors
pred_mask <- app(preds, sum, na.rm = FALSE) # to get NA if any layer has NA
totcl <- mask(totcl, pred_mask)
rm(pred_mask)

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
# the model will be tested on widespread crops 
# the idea is to select crops with the most homogeneous distribution across reg
# defined as a CV < 100% and the minimum abundance across regions must be > 0.05% of total cropland

sc <- fread("AuxData/CropAbundanceSource.csv")
crops <- sc[,ifelse(source == "SPAM", SPAM_Name, Monf_Name)]
rabun <- rast(file.path(indir,sc$file))
names(rabun) <- crops

# Extract crop abundances
d[, (crops):= extract(rabun, cell)]
# set NA abundances to 0
for(j in crops) set(d, i = d[,.I[is.na(.SD)], .SDcols = j], j = j, 0)
# get total abundance for each crop at each region
regtots <- d[, lapply(.SD, sum, na.rm = T), .SDcols = (crops), by = reg]
regtots <- melt(regtots, measure.vars = crops, variable.name = "crop", 
                value.name = "area") # area is in hectares
# crop area fraction of total cropland of each region
regtots[, frac:= area/sum(area), by = reg]

# coefficient of variation of abundance (fraction) across regions
regvar <- regtots[, .(CV = sd(frac)/mean(frac),
                      mfrac = min(frac)),
                  by = crop]
# select testing crops
crops <- regvar[CV < 1 & mfrac > .0005, crop] 
crops <- as.character(crops)

# remove non-selected crops from d and rabun
ctorm <- names(rabun)[!names(rabun) %in% crops]
d[, (ctorm):= NULL]
rabun <- rabun[crops]

# crops by source
monf <- crops[crops %in% sc[source == 'Monf', crop]]
spam <- crops[crops %in% sc[source == 'SPAM', crop]]

# change abundances to fraction
d[, (crops):= lapply(.SD, function(x) x/tot), .SDcols = (crops)]

## Extract predictor values --------
preds_names <- names(preds)
d[, (preds_names):= extract(preds, cell)]


## get DQI  ------------
# it only affects Monfreda Crops. 
dqi.fn <- paste0("Monfreda/DataQualityIndex/", monf, "_DQI.tif")
rdqi <- rast(file.path(indir,dqi.fn))
monf.dqi <- paste0(monf, "_dqi")

## Extract DQI values --------
d[, (monf.dqi):= extract(rdqi, cell)]

# DQI is a measure of environmental variability. Must be changed to an index
# NA means that it came from county level (maximum data quality index)
envCV_to_dqi <- function(x) ifelse(is.na(x), 1, pmax(1 - x/4,0))
for(j in monf.dqi) set(d, j = j, value = envCV_to_dqi(d[[j]]))

# output dir for CV results
tundir <- "OutData/SDM/CV/Tuning"
dir.create(tundir, F, T)


## Training Subset -----------
# 2% of the data will be used for model training. 
# same cells for all crops and models. Sampled regularly
ss = 0.02
tr <- which(d$cell%%(1/ss) == 0)


# TUNE MODELS #################################################################
## 0. NULL MODEL ------
# average crop fraction 
null_rmse <- function(d, crops){
  rmse <- data.table(crops = crops)
  rmse[, null:= mapply(frmse, 
                       d[,.SD, .SDcols = (crops)], 
                       d[,lapply(.SD,mean), .SDcols = (crops)])]
}

null_brmse <- function(d, crops){
  brmse <- data.table(crops = crops)
  brmse[, null:= mapply(fbrmse, 
                        d[,.SD, .SDcols = (crops)], 
                        d[,lapply(.SD,mean), .SDcols = (crops)])]
}


## 1. MAXENT ------------
### 1.1 Method for presence ------
rmse <- null_rmse(d, crops)
brmse <- null_brmse(d, crops)
cors <- data.table(crops = crops)

# Set n# of presence cells 
npr = 2000

set.seed(0)
for(i in 1:length(crops)){
  js <- if(crops[i] %in% monf) 1:7 else c(1:3, 5L, 6L)
  for(j in js){
    for(k in 1:nk){
      pr <- fpr(d[tr], crops[i], j, npr, colk = 'reg', k = k)
      m <- maxnet(c(rep(1, npr), rep(0, d[tr][reg != k, .N])),
                  d[tr][c(pr, which(reg != k)), ..preds_names])  
      lp <- predict(m, d[reg == k, ..preds_names], type = "link")
      mean_abun <- mean(d[[crops[i]]][d$reg == k])
      clog <- clogtr(mean_abun, lp)
      d[reg == k, paste0(crops[i],".ME"):= clog]
      rm(pr, m, lp, clog, mean_abun); gc()
      print(paste0("ME: ", i, ".", j, ".",  k))
    }
    set(rmse, i, paste0("ME.PrMeth.", j),
        frmse(d[[crops[i]]], d[[paste0(crops[i],".ME")]]))
    set(brmse, i, paste0("ME.PrMeth.", j),
        fbrmse(d[[crops[i]]], d[[paste0(crops[i],".ME")]]))
    set(cors, i, paste0("ME.PrMeth.", j),
        cor(d[[crops[i]]], d[[paste0(crops[i],".ME")]]))
  }
  d[, paste0(crops[i], ".ME"):= NULL] #clean up predictions
}

fwrite(rmse, file.path(tundir, "rmse_ME_PresenceMethods.csv"))
fwrite(brmse, file.path(tundir, "brmse_ME_PresenceMethods.csv"))
fwrite(cors, file.path(tundir, "cors_ME_PresenceMethods.csv"))

### 1.2 Shrink -------
rmse <- null_rmse(d, crops)
brmse <- null_brmse(d, crops)
cors <- data.table(crops = crops)

npr = 2000
#j <- which.max(colMeans(cors[,-1]))
j <- 2L
shr <- c(0.5,1,2,3,4,8,12)

set.seed(0)
for(i in 1:length(crops)){
  for(s in 1:length(shr)){
    for(k in 1:nk){
      pr <- fpr(d[tr], crops[i], j, npr, colk = 'reg', k = k)
      m <- maxnet(c(rep(1, npr), rep(0, d[tr][reg != k, .N])),
                  d[tr][c(pr, which(reg != k)), ..preds_names],
                  regmult = shr[s])  
      lp <- predict(m, d[reg == k, ..preds_names], type = "link")
      mean_abun <- mean(d[[crops[i]]][d$reg == k])
      clog <- clogtr(mean_abun, lp)
      d[reg == k, paste0(crops[i],".ME"):= clog]
      rm(pr, m, lp, clog, mean_abun); gc()
      print(paste0("ME: ", i, ".", s, ".",  k))
    }
    set(rmse, i, paste0("ME.regmult.", shr[s]),
        frmse(d[[crops[i]]], d[[paste0(crops[i],".ME")]]))
    set(brmse, i, paste0("ME.regmul.", shr[s]),
        fbrmse(d[[crops[i]]], d[[paste0(crops[i],".ME")]]))
    set(cors, i, paste0("ME.regmult.", shr[s]),
        cor(d[[crops[i]]], d[[paste0(crops[i],".ME")]]))
  }
  d[, paste0(crops[i], ".ME"):= NULL] #clean up predictions
}

fwrite(rmse, file.path(tundir, "rmse_ME_regmult.csv"))
fwrite(brmse, file.path(tundir, "brmse_ME_regmult.csv"))
fwrite(cors, file.path(tundir, "cors_ME_regmult.csv"))

### 1.3 npr -----------
rmse <- null_rmse(d, crops)
brmse <- null_brmse(d, crops)
cors <- data.table(crops = crops)

# for some reason, model performs worse when considering all data
j_pm <- 2L
regmult <- 3.5

# Set n# of presence cells 
npr = seq(2000,12000,2000)
nbg = npr * 5

set.seed(0)
for(i in 1:length(crops)){
  for(n in 1:length(npr)){
    for(k in 1:nk){
      pr <- fpr(d, crops[i], j_pm, npr[n], colk = 'reg', k = k)
      bg <- sample(d[reg != k, .I], nbg[n])
      m <- maxnet(c(rep(1, npr[n]), rep(0, nbg[n])),
                  d[c(pr, bg), ..preds_names],
                  regmult = regmult)  
      lp <- predict(m, d[reg == k, ..preds_names], type = "link")
      mean_abun <- mean(d[[crops[i]]][d$reg == k])
      clog <- clogtr(mean_abun, lp)
      d[reg == k, paste0(crops[i],".ME"):= clog]
      rm(pr, m, lp, clog, mean_abun); gc()
      print(paste0("ME: ", i, ".", n, ".",  k))
    }
    set(rmse, i, paste0("ME.npr.", npr[n]),
        frmse(d[[crops[i]]], d[[paste0(crops[i],".ME")]]))
    set(brmse, i, paste0("ME.npr.", npr[n]),
        fbrmse(d[[crops[i]]], d[[paste0(crops[i],".ME")]]))
    set(cors, i, paste0("ME.npr.", npr[n]),
        cor(d[[crops[i]]], d[[paste0(crops[i],".ME")]]))
  }
}

fwrite(rmse, file.path(tundir, "rmse_ME_npr.csv"))
fwrite(brmse, file.path(tundir, "brmse_ME_npr.csv"))
fwrite(cors, file.path(tundir, "cors_ME_npr.csv"))

### 1.4 npr after subsample -----------
rmse <- null_rmse(d, crops)
brmse <- null_brmse(d, crops)
cors <- data.table(crops = crops)

# combining subsampling regularly before sampling presence observations
j_pm <- 2
regmult <- 3.5


# Set n# of presence cells 
# 2% of the data will be used for model training. 
# same cells for all crops and models. Sampled regularly
snpr <- expand.grid(ss = c(0.01, 0.04, 0.1),
                    npr = c(1000, 4000, 10000))

set.seed(0)
for(i in 1:length(crops)){
  for(n in 1:nrow(snpr)){
    ss <- snpr$ss[n]
    npr <- snpr$npr[n]
    for(k in 1:nk){
      ntr <- which(d$cell%%(1/ss) == 0)
      pr <- fpr(d[ntr], crops[i], j_pm, npr, colk = 'reg', k = k)
      m <- maxnet(c(rep(1, npr), rep(0, d[ntr][reg != k, .N])),
                  d[ntr][c(pr, which(reg != k)), ..preds_names],
                  regmult = regmult)  
      lp <- predict(m, d[reg == k, ..preds_names], type = "link")
      mean_abun <- mean(d[[crops[i]]][d$reg == k])
      clog <- clogtr(mean_abun, lp)
      d[reg == k, paste0(crops[i],".ME"):= clog]
      rm(pr, m, lp, clog, mean_abun); gc()
      print(paste0("ME: ", i, ".", n, ".",  k))
    }
    set(rmse, i, paste0("ME.snpr.", n),
        frmse(d[[crops[i]]], d[[paste0(crops[i],".ME")]]))
    set(brmse, i, paste0("ME.snpr.", n),
        fbrmse(d[[crops[i]]], d[[paste0(crops[i],".ME")]]))
    set(cors, i, paste0("ME.snpr.", n),
        cor(d[[crops[i]]], d[[paste0(crops[i],".ME")]]))
  }
  d[, paste0(crops[i], ".ME"):= NULL] #clean up predictions
}

fwrite(rmse, file.path(tundir, "rmse_ME_snr.csv"))
fwrite(brmse, file.path(tundir, "brmse_ME_snpr.csv"))
fwrite(cors, file.path(tundir, "cors_ME_snpr.csv"))

# the smallest subsample size and number of presence returns the smallest rmse
# but min bias rmse and max cor is with ss = 0.04 and npr = 4000


## 2. RF ---------------
rmse <- null_rmse(d, crops)
brmse <- null_brmse(d, crops)
cors <- data.table(crops = crops)

tp <- expand.grid(
  mtry = c(1:4,seq(6,nlyr(preds),4)),
  min.node.size = seq(3,15,3)
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
    set(brmse, i, paste0("RF.tp.", j),
        fbrmse(d[[crops[i]]], d[[paste0(crops[i],".RF")]]))
    set(cors, i, paste0("RF.tp.", j),
        cor(d[[crops[i]]], d[[paste0(crops[i],".RF")]]))
  }
  d[, paste0(crops[i], ".RF"):= NULL] #clean up predictions
}

fwrite(rmse, file.path(tundir, "rmse_RF_tun.csv"))
fwrite(brmse, file.path(tundir, "brmse_RF_tun.csv"))
fwrite(cors, file.path(tundir, "cors_RF_tun.csv"))
fwrite(tp, file.path(tundir, "RF_tun_tp.csv"))

# max cor and min rmse is obtained with mtry = 1 and min.node.size = 12


## 3. BRT ---------------

### search 1 ----
rmse <- null_rmse(d, crops)
brmse <- null_brmse(d, crops)
cors <- data.table(crops = crops)

tp <- expand.grid(
  eta = c(.01, .1),
  max_depth = c(3, 5, 7),
  min_child_weight = c(1, 5, 10),
  subsample = c(.65, .8, 1), 
  colsample_bytree = 1,
  nrounds = seq(1200,1800,300)              
)

set.seed(0)
for(i in 1:length(crops)){  
  for(j in 1:nrow(tp)){
    for(k in 1:nk){
      params <- list(eta = tp$eta[j],
                     max_depth = tp$max_depth[j],
                     min_child_weight = tp$min_child_weight[j],
                     subsample = tp$subsample[j],
                     colsample_bytree = tp$colsample_bytree[j]
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
    set(brmse, i, paste0("BRT.tp.", j),
        fbrmse(d[[crops[i]]], d[[paste0(crops[i],".BRT")]]))
    set(cors, i, paste0("BRT.tp.", j),
        cor(d[[crops[i]]], d[[paste0(crops[i],".BRT")]]))
  }
  d[, paste0(crops[i], ".BRT"):= NULL] #clean up predictions
}

fwrite(rmse, file.path(tundir, "rmse_BRT_tun.csv"))
fwrite(brmse, file.path(tundir, "brmse_BRT_tun.csv"))
fwrite(cors, file.path(tundir, "cors_BRT_tun.csv"))
fwrite(tp, file.path(tundir, "BRT_tun_tp.csv"))

### search 2 ---------
# lower eta 
# lower max_depth
# higher min_child_weight
# subsample effect unclear
# less rounds
# col sample by tree?
rmse <- null_rmse(d, crops)
brmse <- null_brmse(d, crops)
cors <- data.table(crops = crops)

tp2 <- expand.grid(
  eta = c(.005, .01, 0.02),
  max_depth = 1:3,
  min_child_weight = c(10,15),
  subsample = 1, 
  colsample_bytree = c(0.8,1),
  nrounds = c(1000,1200)              
)

set.seed(0)
for(i in 1:length(crops)){  
  for(j in 1:nrow(tp2)){
    for(k in 1:nk){
      params <- list(eta = tp2$eta[j],
                     max_depth = tp2$max_depth[j],
                     min_child_weight = tp2$min_child_weight[j],
                     subsample = tp2$subsample[j],
                     colsample_bytree = tp2$colsample_bytree[j]
      )
      brt <- xgboost(params = params,
                     data = as.matrix(d[tr][reg != k, ..preds_names]), 
                     label = d[tr][reg != k][[crops[i]]],
                     nrounds = tp2$nrounds[j], 
                     verbose = 0)
      pv <- predict(brt, as.matrix(d[reg == k, ..preds_names]))
      d[reg == k, paste0(crops[i],".BRT"):= pv]
      print(paste0("BRT: ", i, ".", j, ".",  k))
    }
    set(rmse, i, paste0("BRT.tp2.", j),
        frmse(d[[crops[i]]], d[[paste0(crops[i],".BRT")]]))
    set(brmse, i, paste0("BRT.tp2.", j),
        fbrmse(d[[crops[i]]], d[[paste0(crops[i],".BRT")]]))
    set(cors, i, paste0("BRT.tp2.", j),
        cor(d[[crops[i]]], d[[paste0(crops[i],".BRT")]]))
  }
  d[, paste0(crops[i], ".BRT"):= NULL] #clean up predictions
}

fwrite(rmse, file.path(tundir, "rmse_BRT_tun2.csv"))
fwrite(brmse, file.path(tundir, "brmse_BRT_tun2.csv"))
fwrite(cors, file.path(tundir, "cors_BRT_tun2.csv"))
fwrite(tp2, file.path(tundir, "BRT_tun_tp2.csv"))

# max_depth set to 2
# explore larger min_child weight (20,40, 60)
# explore colsample_bytree (0.4, 0.6, 0.8)

# explore slightly larger eta (>0.02 but <0.1) (seq(0.02,0.08,0.03))
# try more nrounds (800, 1000, 1200, 1400)

### search 3 ---------
rmse <- null_rmse(d, crops)
brmse <- null_brmse(d, crops)
cors <- data.table(crops = crops)

tp3 <- expand.grid(
  eta = seq(.02, .08, 0.03),
  max_depth = 2,
  min_child_weight = c(20,40,60),
  subsample = 1, 
  colsample_bytree = c(.4,.6,.8),
  nrounds = seq(800,1400,200)              
)

set.seed(0)
for(i in 1:length(crops)){  
  for(j in 1:nrow(tp3)){
    for(k in 1:nk){
      params <- list(eta = tp3$eta[j],
                     max_depth = tp3$max_depth[j],
                     min_child_weight = tp3$min_child_weight[j],
                     subsample = tp3$subsample[j],
                     colsample_bytree = tp3$colsample_bytree[j]
      )
      brt <- xgboost(params = params,
                     data = as.matrix(d[tr][reg != k, ..preds_names]), 
                     label = d[tr][reg != k][[crops[i]]],
                     nrounds = tp3$nrounds[j], 
                     verbose = 0)
      pv <- predict(brt, as.matrix(d[reg == k, ..preds_names]))
      d[reg == k, paste0(crops[i],".BRT"):= pv]
      print(paste0("BRT: ", i, ".", j, ".",  k))
    }
    set(rmse, i, paste0("BRT.tp3.", j),
        frmse(d[[crops[i]]], d[[paste0(crops[i],".BRT")]]))
    set(brmse, i, paste0("BRT.tp3.", j),
        fbrmse(d[[crops[i]]], d[[paste0(crops[i],".BRT")]]))
    set(cors, i, paste0("BRT.tp3.", j),
        cor(d[[crops[i]]], d[[paste0(crops[i],".BRT")]]))
  }
  d[, paste0(crops[i], ".BRT"):= NULL] #clean up predictions
}

fwrite(rmse, file.path(tundir, "rmse_BRT_tun3.csv"))
fwrite(brmse, file.path(tundir, "brmse_BRT_tun3.csv"))
fwrite(cors, file.path(tundir, "cors_BRT_tun3.csv"))
fwrite(tp3, file.path(tundir, "BRT_tun_tp3.csv"))

# set max_depth to 2
# set eta to 0.02
# explore smaller colsample_bytree (0.11, 0.22, 0.4)
# explore smaller nrounds
# min_child_weight doesn't have a strong effect, set to 50

### search 4 ---------
rmse <- null_rmse(d, crops)
brmse <- null_brmse(d, crops)
cors <- data.table(crops = crops)

tp4 <- expand.grid(
  eta = .02,
  max_depth = 2,
  min_child_weight = 50,
  subsample = 1, 
  colsample_bytree = c(.11,.22,.4),
  nrounds = seq(400,800,200)              
)

set.seed(0)
for(i in 1:length(crops)){  
  for(j in 1:nrow(tp4)){
    for(k in 1:nk){
      params <- list(eta = tp4$eta[j],
                     max_depth = tp4$max_depth[j],
                     min_child_weight = tp4$min_child_weight[j],
                     subsample = tp4$subsample[j],
                     colsample_bytree = tp4$colsample_bytree[j]
      )
      brt <- xgboost(params = params,
                     data = as.matrix(d[tr][reg != k, ..preds_names]), 
                     label = d[tr][reg != k][[crops[i]]],
                     nrounds = tp4$nrounds[j], 
                     verbose = 0)
      pv <- predict(brt, as.matrix(d[reg == k, ..preds_names]))
      d[reg == k, paste0(crops[i],".BRT"):= pv]
      print(paste0("BRT: ", i, ".", j, ".",  k))
    }
    set(rmse, i, paste0("BRT.tp4.", j),
        frmse(d[[crops[i]]], d[[paste0(crops[i],".BRT")]]))
    set(brmse, i, paste0("BRT.tp4.", j),
        fbrmse(d[[crops[i]]], d[[paste0(crops[i],".BRT")]]))
    set(cors, i, paste0("BRT.tp4.", j),
        cor(d[[crops[i]]], d[[paste0(crops[i],".BRT")]]))
  }
  d[, paste0(crops[i], ".BRT"):= NULL] #clean up predictions
}

fwrite(rmse, file.path(tundir, "rmse_BRT_tun4.csv"))
fwrite(brmse, file.path(tundir, "brmse_BRT_tun4.csv"))
fwrite(cors, file.path(tundir, "cors_BRT_tun4.csv"))
fwrite(tp4, file.path(tundir, "BRT_tun_tp4.csv"))

# max cor is obtained with eta = 0.02, max_depth = 3, min_child_weight = 50, colsample_bytree = 0.4, nrounds = 800, subsample = 1. 



# FINAL MODELS  ############################################################

# maxent parameters
j_pm = 2 ; regmult = 3.5 ; ss = 0.04 ; npr = 4000

# random forest parameters
mtry = 1 ; min.node.size = 12

# boosted regression trees parameters
brt.pars <- list(eta = 0.02,
                 max_depth = 3,
                 min_child_weight = 50,
                 subsample = 1,
                 colsample_bytree = 0.4)
nrounds = 800

set.seed(1)
for(i in 1:length(crops)){
  for(k in 1:nk){
    # maxent
    ntr <- which(d$cell%%(1/ss) == 0)
    pr <- fpr(d[ntr], crops[i], j_pm, npr, colk = 'reg', k = k)
    m <- maxnet(c(rep(1, npr), rep(0, d[ntr][reg != k, .N])),
                d[ntr][c(pr, which(reg != k)), ..preds_names],
                regmult = regmult)  
    lp <- predict(m, d[reg == k, ..preds_names], type = "link")
    mean_abun <- mean(d[[crops[i]]][d$reg == k])
    clog <- clogtr(mean_abun, lp)
    d[reg == k, paste0(crops[i],".ME"):= clog]
    rm(ntr, pr, m, lp, clog, mean_abun); gc()
    print(paste0(i, ".",k, ".ME"))
    
    # Random Forest
    rf <- ranger(x = d[reg != k, ..preds_names], 
                 y = d[reg != k][[crops[i]]],
                 num.trees = 750,
                 mtry = mtry,
                 min.node.size = min.node.size)
    pv <- predict(rf, d[reg == k, ..preds_names])$predictions
    d[reg == k, paste0(crops[i],".RF"):= pv]
    rm(rf, pv); gc()
    print(paste0(i, ".",k, ".RF"))
    
    # Boosted Regression Trees
    brt <- xgboost(params = brt.pars,
                   data = as.matrix(d[reg != k, ..preds_names]), 
                   label = d[reg != k][[crops[i]]],
                   nrounds = nrounds, 
                   verbose = 0)
    pv <- predict(brt, as.matrix(d[reg == k, ..preds_names]))
    d[reg == k, paste0(crops[i],".BRT"):= pv]
    rm(brt, pv); gc()
    print(paste0(i, ".",k, ".BRT"))
  }
}

### Write on Disk ----
mods <- c(".ME", ".RF", ".BRT")
predir <- "OutData/SDM/CV/Model_preds"
dir.create(predir, F, T)

for(i in 1:length(crops)){
  for(j in 1:length(mods)){
    r <- rast(totcl)
    v <- rep(NA_real_, ncell(r))
    ij <- paste0(crops[i], mods[j])
    v[d$cell] <- d[[ij]]
    values(r) <- v
    fn <- file.path(predir, paste0(ij, ".tif"))
    wopts <- list(names = ij, filetype = "GTiff",
                  gdal=c("COMPRESS=Deflate","PREDICTOR=1", "ZLEVEL=6"))
    writeRaster(r, filename = fn, overwrite = T, wopt = wopts)
  }
}


# COMPARE MODELS ###########################################################
### Read from Disk ----
mods <- c(".ME", ".RF", ".BRT")
predir <- "OutData/SDM/CV/Model_preds"
outfn <- Sys.glob(file.path(predir, "*.tif"))
outr <- rast(outfn)

cn <- gsub(paste0(predir, "/"), "", gsub(".tif", "", outfn))
names(outr) <- cn

d[, (preds_names):= NULL]
d[, (monf.dqi):= NULL]
d[, (cn):= extract(outr, cell)]

rmse <- null_rmse(d, crops)
brmse <- null_brmse(d, crops)
cors <- data.table(crops = crops)

for(i in 1:length(crops)){
  for (m in mods){
    set(rmse, i, m, frmse(d[[crops[i]]], d[[paste0(crops[i], m)]]))
    set(brmse, i, m, fbrmse(d[[crops[i]]], d[[paste0(crops[i], m)]]))
    set(cors, i, m, cor(d[[crops[i]]], d[[paste0(crops[i], m)]]))
  }
}

outdir <- "OutData/SDM/CV"
fwrite(cors, file.path(outdir, "Tunned_Model_Comp_cor"))
fwrite(rmse, file.path(outdir, "Tunned_Model_Comp_rmse"))
fwrite(brmse, file.path(outdir, "Tunned_Model_Comp_brmse"))

colMeans(cors[,-1])
colMeans(rmse[,-1])
colMeans(brmse[,-1])

plot(c(rabun$wheat/totcl,
       outr$wheat.BRT,
       outr$wheat.ME, 
       outr$wheat.RF), maxcell = ncell(totcl), mar = c(1,2,1,5))
