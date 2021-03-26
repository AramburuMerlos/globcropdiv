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

# create data.table ------------
d <- data.table(cell = 1:ncell(totcl), tot = values(totcl)[,1])

# remove NA (non cropland) cells, order by cell number
d <- d[!is.na(tot),]


## Select Crop Species for testing------
crops <- c("onion", "tomato", "wheat", "maize")

# Extract crop abundances
sc <- fread("AuxData/CropAbundanceSource.csv")
fn <- sc[data.table(crop = crops), on = .NATURAL]$file
rabun <- rast(file.path(indir,fn))
names(rabun) <- crops

d[, (crops):= extract(rabun, cell)]

# set NA abundances to 0
for(j in crops) set(d, i = d[,.I[is.na(.SD)], .SDcols = j], j = j, 0)

# change abundances to fraction
d[, (crops):= lapply(.SD, function(x) x/tot), .SDcols = (crops)]

## Extract predictor values --------
preds_names <- names(preds)
d[, (preds_names):= extract(preds, cell)]


# TUNE MODELS #################################################################
## Training Subset -----------
ss = 0.02
tr <- which(d$cell%%(1/ss) == 0)

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

## 2. RF ---------------
oobe <- data.table(crops = crops)

tp <- expand.grid(
  mtry = c(1,2,4),
  min.node.size = seq(20,100,10),
  sample.fraction = c(.25, .5, .75)
)

set.seed(0)
for(i in 1:length(crops)){  
  for(j in 1:nrow(tp)){
    rf <- ranger(x = d[tr, ..preds_names], 
                 y = d[[crops[i]]],
                 num.trees = 500,
                 mtry = tp$mtry[j],
                 min.node.size = tp$min.node.size[j],
                 sample.fraction = tp$sample.fraction[j])
    set(oobe, i, paste0("RF.tp.", j), sqrt(rf$prediction.error))
    print(paste0("RF: ", i, ".", j))
  }
}

apply(oobe[,-1], 1, which.min)
tp[26,]
tp[23,]
tp[22,]
which.min(colMeans(oobe[,-1]))
rf.tp <- tp[which.min(colMeans(oobe[,-1])),]

## 3. BRT ---------------
min_rmse <- data.table(crops = crops)
opt_tree <- data.table(crops = crops)

tp1 <- expand.grid(
  eta = c(.01, .1),
  max_depth = c(1, 2, 4),
  min_child_weight = c(1, 10, 20),
  subsample = c(.2, .5, 1), 
  colsample_bytree = c(.11, 0.5, 1)
)

tp2 <- expand.grid(
  eta = c(0.05,.01, 0.2, 0.3),
  max_depth = c(4, 8, 16, 32),
  min_child_weight = c(1),
  subsample = c(0.5), 
  colsample_bytree = c(1)
)

set.seed(0)
for(i in 1:length(crops)){  
  for(j in 1:nrow(tp2)){
    params <- list(eta = tp2$eta[j],
                   max_depth = tp2$max_depth[j],
                   min_child_weight = tp2$min_child_weight[j],
                   subsample = tp2$subsample[j],
                   colsample_bytree = tp2$colsample_bytree[j]
    )
    brt <- xgb.cv(params = params,
                  data = as.matrix(d[tr, ..preds_names]), 
                  label = d[tr][[crops[i]]],
                  nrounds = 2000,
                  nfold = 5,
                  verbose = 0, 
                  early_stopping_rounds = 10)
    set(opt_tree, i, paste0("BRT.tp.", j), 
        which.min(brt$evaluation_log$test_rmse_mean))
    set(min_rmse, i, paste0("BRT.tp.", j), 
        min(brt$evaluation_log$test_rmse_mean))
    print(paste0("BRT: ", i, ".", j))
  }
}

apply(min_rmse[,-1], 1, which.min)
which.min(colMeans(min_rmse[,-1]))
opt_tree
min_rmse



# FINAL MODELS  ############################################################

# maxent parameters from regional CV
j_pm = 2 ; regmult = 3.5 ; ss = 0.04 ; npr = 4000

# random forest parameters
mtry = 1 ; min.node.size = 90; sample.fraction = 0.25

# boosted regression trees parameters
nrounds = 650
brt.pars <- list(eta = 0.01,
                 max_depth = 32,
                 min_child_weight = 1,
                 subsample = .5,
                 colsample_bytree = 1)

set.seed(1)

# k folds
nk <- 5
ks <- c(rep(1:nk, floor(d[,.N]/nk)), 1:(d[,.N] %% nk))
d[, fold:= sample(ks)]

for(i in 1:length(crops)){
  for(k in 1:nk){
    # maxent
    ntr <- which(d$cell%%(1/ss) == 0)
    pr <- fpr(d[ntr], crops[i], j_pm, npr, colk = 'fold', k = k)
    m <- maxnet(c(rep(1, npr), rep(0, d[ntr][fold != k, .N])),
                d[ntr][c(pr, which(fold != k)), ..preds_names],
                regmult = regmult)  
    lp <- predict(m, d[fold == k, ..preds_names], type = "link")
    mean_abun <- mean(d[[crops[i]]][d$fold == k])
    clog <- clogtr(mean_abun, lp)
    d[fold == k, paste0(crops[i],".ME"):= clog]
    rm(ntr, pr, m, lp, clog, mean_abun); gc()
    print(paste0(i, ".", k, ".ME"))
    
    # Random Forest
    rf <- ranger(x = d[fold != k, ..preds_names], 
                 y = d[fold != k][[crops[i]]],
                 num.trees = 1000,
                 mtry = mtry,
                 min.node.size = min.node.size, 
                 sample.fraction = sample.fraction)
    pv <- predict(rf, d[fold == k, ..preds_names])$predictions
    d[fold == k, paste0(crops[i],".RF"):= pv]
    rm(rf, pv); gc()
    print(paste0(i, ".", k, ".RF"))
    
    # Boosted Regression Trees
    brt <- xgboost(params = brt.pars,
                   data = as.matrix(d[fold != k, ..preds_names]), 
                   label = d[fold != k][[crops[i]]],
                   nrounds = nrounds)
    pv <- predict(brt, as.matrix(d[fold == k, ..preds_names]))
    d[fold == k, paste0(crops[i],".BRT"):= pv]
    rm(brt, pv); gc()
    print(paste0(i, ".", k, ".BRT"))
  }
}

### Write on Disk ----
mods <- c(".ME", ".RF", ".BRT")
predir <- "OutData/SDM/NormalCV/Model_preds"
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
predir <- "OutData/SDM/NormalCV/Model_preds"
outfn <- Sys.glob(file.path(predir, "*.tif"))
outr <- rast(outfn)

cn <- gsub(paste0(predir, "/"), "", gsub(".tif", "", outfn))
names(outr) <- cn

d[, (preds_names):= NULL]
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

outdir <- "OutData/SDM/NormalCV"
fwrite(cors, file.path(outdir, "Tunned_Model_Comp_cor.csv"))
fwrite(rmse, file.path(outdir, "Tunned_Model_Comp_rmse.csv"))
fwrite(brmse, file.path(outdir, "Tunned_Model_Comp_brmse.csv"))

colMeans(cors[,-1])
colMeans(rmse[,-1])
colMeans(brmse[,-1])

par(mfrow = c(2,2))
plot(c(rabun$wheat/totcl,
        outr$wheat.BRT,
        outr$wheat.ME, 
        outr$wheat.RF), 
     maxcell = ncell(totcl), plg = list(shrink = 0.7, cex = 0.7), 
     cex.axis = 0.8, mar = c(1,2,1,6), axes = FALSE, range = c(0,1))
