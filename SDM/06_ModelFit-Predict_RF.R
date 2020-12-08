library(raster)
library(terra)
library(data.table)
library(randomForest)

nwd <- nchar(getwd())
if(substr(getwd(),  nwd - 10, nwd) != "globcropdiv") warning("See 0000_wd.R")

if(dir.exists("D:/globcropdiv")){
  wcpath <- "D:/WorldClim/2.1/wc5min"
  sgpath <- "D:/SoilGrids"
  gdpath <- "D:/globcropdiv"
  aqpath <- "D:/AQUASTAT"
  monpath <- "D:/Monfreda"
} else {
  wcpath <- "InData/WorldClim/2.1/wc5min"
  sgpath <- "InData/SoilGrids"
  aqpath <- "InData/AQUASTAT"
  gdpath <- "OutData"
  monpath <- "InData/Monfreda"
}

# Predictors --------------
preds <- c(rast(Sys.glob(file.path(wcpath, "bioc/*.tif"))), # WorldClim Bioclim
           rast(file.path(wcpath, "extra/GDD.tif")), 
           rast(file.path(wcpath, "extra/AI.tif")),  
           rast(file.path(wcpath, "extra/PET_seasonality.tif")),  
           rast(file.path(sgpath, "phh2o/phh2o_0-15cm_mean_5min.tif")), 
           rast(file.path(aqpath, "gmia_v5_aei_pct.asc"))  # irrigation
)


# Crop Abundance Data ----------
fn <- Sys.glob(file.path(monpath, "GeoTiff/*/*HarvestedAreaFraction.tif"))
crops <- sub(file.path(monpath, "GeoTiff/*"), "", 
             sub("/*_HarvestedAreaFraction.tif", "", fn))
crops <- substr(crops, 1, ceiling(nchar(crops)/2)-1)
abund <- rast(fn)
# compare predictors and abundance Geometries
compareGeom(preds, abund, lyrs = F)

# subsample size --------
# eventually all presence data will be used to train the model(s)
# and only absence (background) data might be subsampled
# (all crops have >> 1e5 absences, maybe keep 2e4 or something like that)

# in the meantime, let's fit the RF with a "small" subsample
maxN <- 2e4 # max number of training obs 
minBg <- 1e3 # min number or background pts 
train_cells <- vector(mode = "list", length = length(crops))
names(train_cells) <- crops

outpath <- file.path(gdpath, "SDM")
dir.create(outpath)
# loop for all crops  ----
outpath <- file.path(outpath, "RandomForest")
dir.create(outpath)

set.seed(1234)
for(i in 1:length(crops)){
  d <- data.table(cell = (1:ncell(abund[[i]])),
                  frac = getValues(raster(abund[[i]])))
  d <- d[!is.na(frac),]
  totPr <- sum(d$frac > 0)# total # of presence values
  if(totPr > (maxN - minBg)){
    d0 <- d[frac == 0,]
    d0 <- d0[sample(1:nrow(d0), minBg),]
    dPr <- d[frac > 0, ]
    dPr[, w:= frac/max(frac)*.9] # *.9 to aviod probabilities of 1
    dPr[, w:= max(w, 0.3)] # set min p to 0.3 (high frac up to x3 more likely)
    dPr <- dPr[sample(1:nrow(dPr), maxN - minBg, prob = w),]
    dPr[, w:=NULL]
  } else if (totPr > maxN/2){
    d0 <- d[frac == 0,]
    d0 <- d0[sample(1:nrow(d0), (maxN - totPr)),] 
    dPr <- d[frac > 0, ]
  } else {
    d0 <- d[frac == 0,]
    d0 <- d0[sample(1:nrow(d0), max(totPr, minBg)),] 
    dPr <- d[frac > 0, ]
  }
  d <- rbind(d0, dPr)
  setorder(d, cell)
  train_cells[[i]] <- d$cell
  d[, cell:= NULL]
  dpreds <- as.data.table(extract(preds, train_cells[[i]]))
  d <- cbind(d, dpreds)
  d <- d[complete.cases(d),]
  rm(dPr, d0, dpreds)
  
  # Random Forest model -------
  # first tune mtry, then run the model. All other parameter are by default
  dtune <- d[sample.int(5000),]
  trf <- tuneRF(dtune[,-"frac"],dtune$frac, stepFactor = 1.5, improve = 0.01,
                trace = FALSE, plot = FALSE)
  mt <- trf[which.min(trf[,2]), 1]
  rm(dtune)
  rfm <- randomForest(frac ~ ., mtry = mt, data = d) 
  saveRDS(rfm, file.path(outpath, paste0("RF", crops[i], ".RDS")))
  
  # Predict ---------
  predict(preds, rfm, na.rm = T, 
          filename = file.path(outpath, paste0(crops[i], "_RFpred.tif")), 
          overwrite = T, 
          wopt = list(names = crops[i], filetype = "GTiff",
                      gdal = c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))
}


r = rast(file.path(outpath, paste0(crops[i], "_RFpred.tif")))
plot(r)
plot(abund[[i]])





