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
               raster(file.path(aqpath, "gmia_v5_aei_pct.asc")))
bioc <- stack(lapply(Sys.glob(file.path(wcpath, "bioc/*.tif")), raster))
preds <- stack(preds, bioc)
cmask <- raster(file.path(gdpath, "Monf_Cropland_Mask.tif"))
preds[[1]] <- mask(preds[[1]], cmask)

# Crop Abundance Data at County Level----------
fn <- Sys.glob(file.path(monpath, "*/*_CountyLevel.tif"))
crops <- sub(paste0(monpath, "/OnlyCounty/"), "", 
             sub("/*_CountyLevel.tif", "", fn))
abund <- stack(lapply(fn, raster))

# compare predictors and abundance Geometries
compareRaster(preds, abund)

# SET THRESHOLD  -----------
# define proportion of abundance data to keep
# i.e. consider as present the th {0,1} most abundant cells
th <- 0.8 # I can play with this value

outpath <- file.path(gdpath, "SDM/Maxent")
dir.create(outpath)
train_cells <- vector(mode = "list", length = length(crops))

maxent()

# loop for all crops  ----
set.seed(1234)
for(i in 1:length(crops)){
  if(maxValue(abund[[i]] > 0)){
    d <- data.table(cell = (1:ncell(abund[[i]])),
                    frac = getValues(abund[[i]]))
    th_frac <- quantile(d$frac, th, na.rm = T, names = F)
    d <- d[frac > th_frac, ]
    pts <- xyFromCell(abund[[i]], d$cell) 
    train_cells[[i]] <- d$cell
    rm(d)
    
    # Maxent Model -------
    m <- maxent(preds, pts, path = file.path(outpath, crops[i]))
    
    # Predict ---------
    predict(m, preds, 
            filename = file.path(outpath, paste0(crops[i], "_Maxent_pred.tif")), 
            overwrite = T, progress = 'text',
            options = c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6"))
  }
}


r = raster(file.path(outpath, paste0(crops[i], "_Maxent_pred.tif")))
plot(r)
plot(abund[[i]])





