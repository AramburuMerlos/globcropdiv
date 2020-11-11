library(data.table)
library(raster)
options(java.parameters = "-Xmx20g" )
library(dismo)

nwd <- nchar(getwd())
if(substr(getwd(),  nwd - 10, nwd) != "globcropdiv") stop("See 0000_wd.R")

if(dir.exists("D:/globcropdiv")){
  wcpath <- "D:/WorldClim/2.1/wc5min"
  sgpath <- "D:/SoilGrids"
  gdpath <- "D:/globcropdiv"
  aqpath <- "D:/AQUASTAT"
  cpath <- "D:/WorldCropAbundance"
} else {
  wcpath <- "InData/WorldClim/2.1/wc5min"
  sgpath <- "InData/SoilGrids"
  aqpath <- "InData/AQUASTAT"
  gdpath <- "OutData"
  cpath <- "InData/WorldCropAbundance"
}

# Predictors --------------
preds <- stack(raster(file.path(sgpath, "phh2o/phh2o_0-15cm_mean_5min.tif")), 
               raster(file.path(wcpath, "extra/GDD.tif")), 
               raster(file.path(wcpath, "extra/AI.tif")),  
               raster(file.path(wcpath, "extra/PET_seasonality.tif")),  
               raster(file.path(aqpath, "gmia_v5_aei_pct.asc")))
bioc <- stack(lapply(Sys.glob(file.path(wcpath, "bioc/*.tif")), raster))
preds <- stack(preds, bioc)
cmask <- raster(file.path(cpath, "Total_Cropland_ha.tif"))
preds[[1]] <- mask(preds[[1]], cmask)

# Crop Abundance Data ----------
fn <- Sys.glob(file.path(cpath, "fraction/*.tif"))
abund <- stack(lapply(fn, raster))
crops <- sub(paste0(cpath, "/fraction/"), "", 
             sub("_...._fr.tif", "", fn))
# compare predictors and abundance Geometries
compareRaster(preds, abund)

# SET "SAMPLING EFFORT"  -----------
# define proportion|number of abundance data (>0) to keep 
# then, sample abundance data using crop fraction as weight 

# sampling effort is set as the minimum proportion that ensures that all (99%) presence cells have a suitability > 0 

outpath <- file.path(gdpath, "SDM/Maxent")
dir.create(outpath)
sampling_effort <- vector(length = length(crops))
train_cells <- vector(mode = "list", length = length(crops))

maxent()

seff <- seq(0.2,1,0.2)

# loop for all crops  ----
set.seed(1234)
for(i in 1:length(crops)){
  for(s in seff){
    d <- data.table(cell = (1:ncell(abund[[i]])),
                    frac = getValues(abund[[i]]))
    d <- d[!is.na(frac),]
    scells <- sample(d$cell, ceiling(nrow(d)*s), prob = d$frac)
    pts <- xyFromCell(abund[[i]], scells) 
    train_cells[[i]] <- scells
    rm(d)
    # fit
    m <- maxent(preds, pts, path = file.path(outpath, crops[i]))
    # Predict
    mp <- predict(m, preds, 
            filename = file.path(outpath, paste0(crops[i], "_Maxent_pred.tif")), 
            overwrite = T, 
            options = c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6"))
    mpm <- mask(mp, abund[[i]])
    mpmv <- getValues(mpm)
    mpmv <- mpmv[!is.na(mpmv)]
    if(sum(mpmv > 0)/ length(mpmv) > 0.99) break
  }
  sampling_effort[i] <- s
  train_cells[[i]] <- scells
}
  
fwrite(data.table(crop = crops, sampling_effort = sampling_effort), 
       file.path(outpath, "sampling_effort.csv"))

saveRDS(train_cells, file.path(outpath, "train_cells.RDS"))



