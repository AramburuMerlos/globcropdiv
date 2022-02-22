library(magrittr)
library(data.table)
library(raster)
library(terra)
library(rnaturalearth)

if(system('hostname', TRUE) %in% c("ESP-RH-9891", "LAPTOP-ST129J47")){
  setwd("D:/globcropdiv/")
} # else if { ... 

sc <- fread("AuxData/CropCategories.csv")

# Data Quality Index (Monf Data) ########################################

# upload Monfreda crop abundance data
mcrops <- sc[source == "Monf", crop] %>%
  paste0(".tif") %>%
  file.path("InData/CropAbundance", .) %>%
  rast()

# total cropland
totcl <- rast("InData/TotalCropland.tif")

# Country and State CV ----
# rast with countries 
spcountries <- ne_countries()
spcountries$cID <- 1:nrow(spcountries)
emptyr <- rast(mcrops[[1]])
rcountries <- rasterize(vect(spcountries), emptyr, field = "cID", touches=TRUE) 

# rast with states / provinces
spstates <- ne_states()
spstates$sID <- 1:nrow(spstates)
rstates <- rasterize(vect(spstates), emptyr, field = "sID", touches=TRUE)

# Growing Degree Days and Precipitation
rain <- rast("InData/WorldClim/2.1/wc5min/bioc/wc2.1_5m_bio_12.tif")
GDD <- rast("InData/WorldClim/2.1/wc5min/extra/GDD.tif")

# mask cropland
rain <- mask(rain, totcl)
GDD <- mask(GDD, totcl)
# for rain, mask irrigated land
# (to avoid computing the CV for highly irrigated areas)
irr <- rast("InData/AQUASTAT/gmia_v5_aei_pct.asc")
irr_m <- arith(irr, function(x) x > 50)
rain <- mask(rain, irr_m, maskvalue = 1)
names(rain) <- "rain"

# create data table with values to compute CV
DT <- data.table(values(c(rcountries, rstates, GDD, rain)))
DT <- DT[complete.cases(DT),]

# add 100mm to all annual precipitation data to avoid high CV because of low Pp
DT[, rain:= rain + 100]

# compute CV by country, get max CV and create rast of climate var per country
cvf <- function(x)sd(x, na.rm = T)/mean(x, na.rm = T)
CVc <- DT[, lapply(.SD, cvf), .SDcols = c("GDD", "rain"), by = "cID"]
CVc[, maxCV := pmax(GDD, rain)]
rcm_country <- as.matrix(CVc[,.(cID,maxCV)])
rcvc <- classify(rcountries, rcm_country, othersNA = TRUE)


# compute CV by state, get max CV and create rast of climate var per state
CVs <- DT[, lapply(.SD, cvf), .SDcols = c("GDD", "rain"), by = "sID"]
CVs[, maxCV := pmax(GDD, rain)]
rcm_state <- as.matrix(CVs[,.(sID, maxCV)])
rcvs <- classify(rstates, rcm_state, othersNA = TRUE)

# DQI rasters -----
# create crop specific climate heterogeneity indexes of input data
# Monfreda data quality
dir1 <- "D:/Monfreda/HarvestedAreaYield175Crops_Geotiff/"
dir2 <- "_HarvAreaYield_Geotiff/"
dir3 <- "_DataQuality_HarvestedArea.tif"

mdq <- sc[source == "Monf", crop] %>%
  paste0(dir1, .,dir2 , ., dir3) %>%
  rast()

# function to select CV depending on Monf DQ
# heterogeneity = 0 for county level data
dqfun <- function(dq, cvc, cvs){
  ifelse(dq == 1, 0, 
         ifelse(dq %in% c(0.75, 0.5), cvs,
                ifelse(dq == 0.25, cvc, NA)))
}

# create DQI raster and save them to disk
outdir <- file.path("InData/CropAbundance/DataQualityIndex")
dir.create(outdir, F, T)

for(i in 1:nlyr(mdq)){
  ci <- names(mcrops)[i]
  dqi <- lapp(c(mdq[[i]], rcvc, rcvs), fun = dqfun,
              filename = file.path(outdir, paste0(ci, "_DQI.tif")), 
              overwrite = T, 
              wopt = list(names = ci, filetype = "GTiff",
                          gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))
}

