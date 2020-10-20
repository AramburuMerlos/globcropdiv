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
  monpath <- "D:/Monfreda"
} else {
  wcpath <- "InData/WorldClim/2.1/wc5min"
  sgpath <- "InData/SoilGrids"
  aqpath <- "InData/AQUASTAT"
  gdpath <- "OutData"
  monpath <- "InData/Monfreda"
}

# Predictors --------------
preds <- stack(raster(file.path(sgpath, "phh2o/phh2o_0-15cm_mean_5min.tif")), 
               raster(file.path(wcpath, "extra/GDD.tif")), 
               raster(file.path(wcpath, "extra/AI.tif")),  
               raster(file.path(wcpath, "extra/PET_seasonality.tif")),  
               raster(file.path(aqpath, "gmia_v5_aei_pct.asc")),
               lapply(Sys.glob(file.path(wcpath, "bioc/*.tif")), raster))
cmask <- raster(file.path(gdpath, "Monf_Cropland_Mask.tif"))
preds[[1]] <- mask(preds[[1]], cmask)

# Crop Abundance Data at County Level----------
fn <- Sys.glob(file.path(monpath, "*/*_CountyLevel.tif"))
crops <- sub(paste0(monpath, "/CountyLevel/"), "", 
             sub("/*_CountyLevel.tif", "", fn))
abund <- stack(lapply(fn, raster))


# compare predictors and abundance Geometries
compareRaster(preds, abund)

# SET THRESHOLD  -----------
# define threshold to identify a crop as present
th <- 0 # I can play with this value

# subsample size --------
# eventually all presence data will be used to train the model(s)

# in the meantime, let's fit Maxent with a "small" subsample
np <- 5e4 # max number of training obs 

outpath <- file.path(gdpath, "Maxent")
dir.create(outpath)
train_cells <- vector(mode = "list", length = length(crops))

maxent()

# loop for all crops  ----
set.seed(1234)
for(i in 1:length(crops)){
  d <- data.table(cell = (1:ncell(abund[[i]])),
                  frac = getValues(abund[[i]]))
  d <- d[frac > th, ]
  if(nrow(d) > np){
    d[, w:= frac/max(frac)*.9] # *.9 to aviod probabilities of 1
    d[, w:= max(w, 0.3)] # set min p to 0.3 (high frac up to x3 more likely)
    d <- d[sample(1:nrow(d), np, prob = w),]
  } 
  pts <- xyFromCell(abund[[i]], d$cell) 
  train_cells[[i]] <- d$cell
  rm(d)
  
  # Maxent Model -------
  m <- 
  saveRDS(m, file.path(outpath, paste0("Maxent_", crops[i], ".RDS")))
  
  # Predict ---------
  predict(preds, m, na.rm = T, 
          filename = file.path(outpath, paste0(crops[i], "_Maxent_pred.tif")), 
          overwrite = T, 
          wopt = list(names = crops[i], filetype = "GTiff",
                      gdal = c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))
}


r = rast(file.path(outpath, paste0(crops[i], "_Maxent_pred.tif")))
plot(r)
plot(abund[[i]])





