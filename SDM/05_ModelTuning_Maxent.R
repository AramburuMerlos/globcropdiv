library(magrittr)
library(data.table)
library(terra)
library(maxnet)

nwd <- nchar(getwd())
if(substr(getwd(),  nwd - 10, nwd) != "globcropdiv") warning("See 0000_wd.R")

# directory of input data
indir <- "D:"

# upload functions
source("Functions/maxent_output_transformation.R")


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


## Set n# of presence and bg cells ------
npr = 2000
nbg = 20000

## get DQI for methods 3 and 4 ------------
# it only affects Monfreda Crops. (SPAM results are the same as in model 2)
dqi.fn <- paste0("Monfreda/DataQualityIndex/", monf, "_DQI.tif")
rdqi <- rast(file.path(indir,dqi.fn))
monf.dqi <- paste0(monf, "_dqi")

## Extract DQI values --------
d[, (monf.dqi):= extract(rdqi, cell)]

# DQI is a measure of environmental variability. Must be changed to an index
# NA means that it came from county level (maximum data quality index)
envCV_to_dqi <- function(x) ifelse(is.na(x), 1, pmax(1 - x/4,0))
for(j in monf.dqi) set(d, j = j, value = envCV_to_dqi(d[[j]]))

## Writing setup ----
wopts <- list(names=crops[i], filetype = "GTiff",
              gdal=c("COMPRESS=Deflate","PREDICTOR=1", "ZLEVEL=6"))
# output dir for rasters with predictions
outdir <- "OutData/SDM/Maxent/CVpreds"
dir.create(outdir, recursive = T)

# RUN MODELS ###############################

#0. NULL MODEL ####
# average crop fraction 
# correlation with null model will be always NA (all same value)
# RMSE with itself is population sd

# 1. BASE MODEL ##########################################################
# Base: data as is. Sample any presence with equal weight
i=k=1
set.seed(01042021)
for(i in 1:length(crops)){
  # regional cross validation
  ldt <- vector("list", nk)
  # sample training data (Method I)
  for(k in 1:nk){ 
    ldt[[k]] <- rbind(
      d[reg != k][sample(.N, nbg)], #background
      d[reg != k & d[[crops[i]]] > 0][sample(.N, npr)] # presence, equal weight
    )
    ldt[[k]][,reg:= k] # ***reg now is the region where it will be tested***
    ldt[[k]][,pa:= rep(c(0L,1L), c(nbg,npr))]
  }
  dtrain <- rbindlist(ldt)
  rm(ldt)
  # fit model
  for(k in 1:nk){
    m <- maxnet(dtrain[reg == k, pa], dtrain[reg == k, ..preds_names])  
    # predicted values (linear predictor)
    lp <- predict(m, d[reg == k, ..preds_names], type = "link")
    # transform linear predictions to cloglog with optimum c 
    mean_abun <- mean(d[[crops[i]]][d$reg == k])
    clog <- clogtr(mean_abun, lp)
    d[reg == k, paste0(crops[i],".m1"):= clog]
    rm(lp, clog, mean_abun)
    print(paste0("m1: crop ", i, "/", length(crops), "; ", k,"/",nk))
  }
}

## Write on Disk ----
for(i in 1:length(crops)){
  r <- rast(totcl)
  v <- rep(NA_real_, ncell(r))
  v[d$cell] <- d[[paste0(crops[i], ".m1")]]
  values(r) <- v
  fn <- file.path(outdir, paste0(crops[i], "_Maxent_CVpred_m1.tif"))
  writeRaster(r, filename = fn, overwrite = T, wopt = wopts)
}



# 2. CROP FRAC as SAMPLE PROB MODEL ############################################# 
# Random sampling of cells where the crop is reported, probability of each cell being sampled is proportional to its crop area (as fraction of cropland area). 
i=k=1
set.seed(01042021)
for(i in 1:length(crops)){
  # regional cross validation
  ldt <- vector("list", nk)
  # sample training data (Method II)
  for(k in 1:nk){ 
    ldt[[k]] <- rbind(
      d[reg != k][sample(.N, nbg)], #background
      d[reg != k][sample(.N, npr, prob = d[[crops[i]]][d$reg != k])]
    )
    ldt[[k]][,reg:= k] # ***reg now is the region where it will be tested***
    ldt[[k]][,pa:= rep(c(0L,1L), c(nbg,npr))]
  }
  dtrain <- rbindlist(ldt)
  rm(ldt)
  # fit model
  for(k in 1:nk){
    m <- maxnet(dtrain[reg == k, pa], dtrain[reg == k, ..preds_names])  
    # predicted values (linear predictor)
    lp <- predict(m, d[reg == k, ..preds_names], type = "link")
    # transform linear predictions to cloglog with optimum c 
    mean_abun <- mean(d[[crops[i]]][d$reg == k])
    clog <- clogtr(mean_abun, lp)
    d[reg == k, paste0(crops[i],".m2"):= clog]
    rm(lp, clog, mean_abun)
    print(paste0("m2: crop ", i, "/", length(crops), "; ", k,"/",nk))
  }
}

## Write on Disk ----
for(i in 1:length(crops)){
  r <- rast(totcl)
  v <- rep(NA_real_, ncell(r))
  v[d$cell] <- d[[paste0(crops[i], ".m2")]]
  values(r) <- v
  fn <- file.path(outdir, paste0(crops[i], "_Maxent_CVpred_m2.tif"))
  writeRaster(r, filename = fn, overwrite = T, wopt = wopts)
}


# 3. DQI MODEL ################################################################
i=k=1
set.seed(01042021)
for(i in 1:length(monf)){
  # regional cross validation
  ldt <- vector("list", nk)
  # sample training data (Method III)
  for(k in 1:nk){ 
    ldt[[k]] <- rbind(
      d[reg != k][sample(.N, nbg)], #background
      d[reg != k][sample(.N, npr, 
                         prob = (d[[monf[i]]] * d[[monf.dqi[i]]])[d$reg != k])]
    )
    ldt[[k]][,reg:= k] # ***reg now is the region where it will be tested***
    ldt[[k]][,pa:= rep(c(0L,1L), c(nbg,npr))]
  }
  dtrain <- rbindlist(ldt)
  rm(ldt)
  # fit model
  for(k in 1:nk){
    m <- maxnet(dtrain[reg == k, pa], dtrain[reg == k, ..preds_names])  
    # predicted values (linear predictor)
    lp <- predict(m, d[reg == k, ..preds_names], type = "link")
    # transform linear predictions to cloglog with optimum c 
    mean_abun <- mean(d[[monf[i]]][d$reg == k])
    clog <- clogtr(mean_abun, lp)
    d[reg == k, paste0(monf[i],".m3"):= clog]
    rm(lp, clog, mean_abun)
    print(paste0("m3: crop ", i, "/", length(monf), "; ", k,"/",nk))
  }
}

## Write on Disk ----
for(i in 1:length(monf)){
  r <- rast(totcl)
  v <- rep(NA_real_, ncell(r))
  v[d$cell] <- d[[paste0(monf[i], ".m3")]]
  values(r) <- v
  fn <- file.path(outdir, paste0(monf[i], "_Maxent_CVpred_m3.tif"))
  writeRaster(r, filename = fn, overwrite = T, wopt = wopts)
}


# 4. BINOM SAMPLING MODEL #####################################################
# presence or absence is determined by a binomial random sampling in which
# the prob of presence is a fun of crop fraction. (p(cell = 1) = crop_fraction)
# the process is repeated until enough presents are detected (might be 1 round)
# here, presence cells will be subsampled to npr
i=k=1
set.seed(01042021)
for(i in 1:length(crops)){
  # regional cross validation
  ldt <- vector("list", nk)
  # define presence with random binomial experiment (Method IV)
  for(k in 1:nk){
    pal <- rep(F, d[reg != k,.N])
    while(sum(pal) < npr){
      pal <- pal | rbinom(d[reg != k,.N], 1, d[[crops[i]]][d$reg != k])  
    }
    ldt[[k]] <- rbind(
      d[reg != k][sample(.N, nbg)], #background
      d[reg != k][pal][sample(.N, npr)] # reduce to npr to compare methods
    )
    ldt[[k]][,reg:= k] # ***reg now is the region where it will be tested***
    ldt[[k]][,pa:= rep(c(0L,1L), c(nbg,npr))]
  }
  # sample training data (Method IV)
  dtrain <- rbindlist(ldt)
  rm(ldt)
  # fit model
  for(k in 1:nk){
    m <- maxnet(dtrain[reg == k, pa], dtrain[reg == k, ..preds_names])  
    # predicted values (linear predictor)
    lp <- predict(m, d[reg == k, ..preds_names], type = "link")
    # transform linear predictions to cloglog with optimum c 
    mean_abun <- mean(d[[crops[i]]][d$reg == k])
    clog <- clogtr(mean_abun, lp)
    d[reg == k, paste0(crops[i],".m4"):= clog]
    rm(lp, clog, mean_abun)
    print(paste0("m4: crop ", i, "/", length(crops), "; ", k,"/",nk))
  }
}

## Write on Disk ----
for(i in 1:length(crops)){
  r <- rast(totcl)
  v <- rep(NA_real_, ncell(r))
  v[d$cell] <- d[[paste0(crops[i], ".m4")]]
  values(r) <- v
  fn <- file.path(outdir, paste0(crops[i], "_Maxent_CVpred_m4.tif"))
  writeRaster(r, filename = fn, overwrite = T, wopt = wopts)
}


#5. BINOM SAMPLING MODEL w DQI ################################################
# same as 4, but the probability of success in the binomial experiment is undermined by the data quality index of that cell
i=k=1
set.seed(01042021)
for(i in 1:length(monf)){
  # regional cross validation
  ldt <- vector("list", nk)
  # define presence with random binomial experiment affected by DQI (Method V)
  for(k in 1:nk){
    pal <- rep(F, d[reg != k,.N])
    while(sum(pal) < npr){
      pal <- pal | rbinom(d[reg != k,.N], 1, 
                          (d[[monf[i]]] * d[[monf.dqi[i]]])[d$reg != k])  
    }
    ldt[[k]] <- rbind(
      d[reg != k][sample(.N, nbg)], #background
      d[reg != k][pal][sample(.N, npr)] # reduce to npr to compare methods
    )
    ldt[[k]][,reg:= k] # ***reg now is the region where it will be tested***
    ldt[[k]][,pa:= rep(c(0L,1L), c(nbg,npr))]
  }
  # sample training data (Method IV)
  dtrain <- rbindlist(ldt)
  rm(ldt)
  # fit model
  for(k in 1:nk){
    m <- maxnet(dtrain[reg == k, pa], dtrain[reg == k, ..preds_names])  
    # predicted values (linear predictor)
    lp <- predict(m, d[reg == k, ..preds_names], type = "link")
    # transform linear predictions to cloglog with optimum c 
    mean_abun <- mean(d[[monf[i]]][d$reg == k])
    clog <- clogtr(mean_abun, lp)
    d[reg == k, paste0(monf[i],".m5"):= clog]
    rm(lp, clog, mean_abun)
    print(paste0("m5: crop ", i, "/", length(monf), "; ", k,"/",nk))
  }
}

## Write on Disk ----
for(i in 1:length(monf)){
  r <- rast(totcl)
  v <- rep(NA_real_, ncell(r))
  v[d$cell] <- d[[paste0(monf[i], ".m5")]]
  values(r) <- v
  fn <- file.path(outdir, paste0(monf[i], "_Maxent_CVpred_m5.tif"))
  writeRaster(r, filename = fn, overwrite = T, wopt = wopts)
}

# READ from Disk #######################################
outfn <- Sys.glob(file.path(outdir, "*.tif"))
outr <- rast(outfn)

cn <- paste0(outdir, "/") %>%
  gsub("", outfn) %>%
  gsub("_Maxent_CVpred_", ".", .) %>%
  gsub(".tif", "", .)
names(outr) <- cn

d[, (cn):= extract(outr, cell)]


# EVAL RESULTS ##########################################
cors <- data.table(crop = crops)

for(i in 1:length(crops)){
  js <- if(crops[i] %in% monf) 1:5 else c(1,2,4)
  cm <- paste0(crops[i], ".m", js)
  cmcors <- sapply(d[, ..cm], function(x)cor(x, d[[crops[i]]], use = "comp"))
  for(j in js) set(cors, i, paste0("m", j), cmcors[paste0(crops[i], ".m", j)])
}
cors
cors[!is.na(m3), lapply(.SD, mean), .SDcols = paste0("m",1:5)]
cors[is.na(m3), lapply(.SD, mean), .SDcols = paste0("m", c(1,2,4))]
fwrite(cors, "Outdata/SDM/Maxent/CVeval_cor.csv")

frmse <- function(x, y) sqrt(mean((x-y)^2, na.rm = T))
sdp <- function(x) sd(x, na.rm = T) * (length(x) - 1)/length(x)
rmse <- data.table(crop = crops)

for(i in 1:length(crops)){
  rmse[crop == crops[i], m0:= sdp(d[[crops[i]]])]
  js <- if(crops[i] %in% monf) 1:5 else c(1,2,4)
  cm <- paste0(crops[i], ".m", js)
  cmrmse <- sapply(d[, ..cm], function(x)frmse(x, d[[crops[i]]]))
  for(j in js) set(rmse, i, paste0("m", j), cmrmse[paste0(crops[i], ".m", j)])
}
rmse
rmse[!is.na(m3), lapply(.SD, mean), .SDcols = paste0("m", 0:5)]
rmse[is.na(m3), lapply(.SD, mean), .SDcols = paste0("m", c(0:2,4))]
fwrite(rmse, "Outdata/SDM/Maxent/CVeval_rmse.csv")
