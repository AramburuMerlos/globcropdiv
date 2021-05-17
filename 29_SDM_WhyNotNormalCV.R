library(magrittr)
library(data.table)
library(terra)
library(xgboost)

if(system('hostname', TRUE) == "ESP-RH-9891"){
  setwd("D:/globcropdiv/")
} else if(system('hostname', TRUE) == "LAPTOP-ST129J47"){ 
  setwd("G:/My Drive/globcropdiv/")
} # else if { ... 

frmse <- function(obs, pred) sqrt(mean((obs-pred)^2, na.rm = T))

# DATA PREP ###########################################################
## Predictors --------------
r_preds <- fread("AuxData/SelectedPreds.csv")$file_name %>% 
  file.path("InData", .) %>% 
  Sys.glob() %>% 
  rast()
preds <- fread("AuxData/SelectedPreds.csv")$short_name

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

## Take sample ------------
# only 5% of the data will be used for this exercise; sampled regularly
d <- d[cell %% (1/.05) == 0]


## Wheat abundance data -----
"InData/CropAbundance/wheat.tif" %>% rast() %>% 
{d[, wheat:= extract(., cell)]}

# change abundances to fraction
d[, wheat:= wheat/tot]

## Extract predictor values --------
d[, (preds):= extract(r_preds, cell)]


# CROSS VALIDATE ######################
# in boosted regression trees, a large number of trees can overfit the data
# that's why it is recommended to select the number of trees with cross validation
# but what if data isn't independent?

## Normal CV ---------------
set.seed(0)

xgb.normCV <- xgb.cv(
  params = list(eta = 0.1, max_depth = 5),
  data = as.matrix(d[, ..preds]),
  label = d[, wheat],
  nrounds = 2000,
  nfold = 5,
  print_every_n = 250,
  objective = "reg:squarederror"
)

## Regional CV ---------------
nrounds <- data.frame(nrounds = seq(0, 2000, 100), test_rmse= NA_real_)

set.seed(0)
for(k in 1:nk){
  dm <- xgb.DMatrix(data = as.matrix(d[reg != k, ..preds]), 
                    label = d[reg != k, wheat])
  mod <- NULL
  for(i in 1:nrow(nrounds)){
    mod <- xgb.train(
      params = list(eta = 0.1, max_depth = 5),
      data = dm,
      nrounds = nrounds$nrounds[i], 
      verbose = 0
    )
    d[reg == k, paste0("y.", i):= predict(mod, as.matrix(.SD)), .SDcols = preds]
    print(paste0(k, "/", nk, " - ", i, "/", nrow(nrounds)))
  }
}

for(i in 1:nrow(nrounds)){
  nrounds$test_rmse[i] <- frmse(d$wheat, d[[paste0("y.", i)]])
}
nrounds

plot(xgb.normCV$evaluation_log$iter,
     xgb.normCV$evaluation_log$test_rmse_mean, 
     type = "line")
lines(nrounds$nrounds, nrounds$test_rmse, col = "red")
which.min(xgb.normCV$evaluation_log$test_rmse_mean)
tail(xgb.normCV$evaluation_log$test_rmse_mean)
