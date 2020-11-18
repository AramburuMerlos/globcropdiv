library(terra)
library(data.table)

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

# Load Predictors and Mask Cropland --------------
preds <- c(rast(Sys.glob(file.path(wcpath, "bioc/*.tif"))), # WorldClim Bioclim
           rast(file.path(wcpath, "extra/GDD.tif")), 
           rast(file.path(wcpath, "extra/AI.tif")),  
           rast(file.path(wcpath, "extra/PET_seasonality.tif")),  
           rast(file.path(sgpath, "phh2o/phh2o_0-15cm_mean_5min.tif")), 
           rast(file.path(aqpath, "gmia_v5_aei_pct.asc"))  # irrigation
)

cmask <- rast(file.path(cpath, "Total_Cropland_ha.tif"))
preds <- mask(preds, cmask)

# Correlation Matrix to identify highly correlated variables
# 
# PCA to see which variables explain most of the first dimensions


