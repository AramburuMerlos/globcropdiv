library(data.table)
library(raster)
library(terra)
library(rnaturalearth)

if(!getwd() %like%  "globcropdiv$") warning("See 0000_wd.R")

sc <- fread("AuxData/CropAbundanceSource.csv")

# Drive or folder where the Raw SPAM and Monfreda data is saved
rd <- "D:/"


# Data Quality Index (Monf Data) ########################################

# upload Monfreda crop abundance data
mfn <- sc[source == "Monf", file]
mcrops <- rast(file.path(rd, mfn))
  
# rast with total cropland (ha) per cell based on Monfreda data
monf_tot <- app(mcrops, fun = sum, na.rm = T)

# raster with 1 for cropland and 0 for non cropland
monf_cropland <- arith(monf_tot, fun = function(x) x > 0)

# Country and State CV ----
# rast with countries 
spcountries <- ne_countries()
spcountries$cID <- 1:nrow(spcountries)
emptyr <- rast(mcrops[[1]])
rcountries <- rasterize(vect(spcountries), emptyr, field = "cID") 

# rast with states / provinces
spstates <- ne_states()
spstates$sID <- 1:nrow(spstates)
rstates <- rasterize(vect(spstates), emptyr, field = "sID")

# Growing Degree Days and Precipitation
rain <- rast(file.path(rd, "WorldClim/2.1/wc5min/bioc/wc2.1_5m_bio_12.tif"))
GDD <- rast(file.path(rd, "WorldClim/2.1/wc5min/extra/GDD.tif"))

# mask cropland
rain <- mask(rain, monf_cropland, maskvalue = 0)
GDD <- mask(GDD, monf_cropland, maskvalue = 0)
# for rain, mask irrigated land
# (to avoid computing the CV for highly irrigated areas)
irr <- rast(paste0(rd, "AQUASTAT/gmia_v5_aei_pct.asc"))
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
mfn_dq <- vector(length = length(mcode))
for(i in 1:length(mcode)){
  mfn_dq[i] <- paste0(rd, "Monfreda/HarvestedAreaYield175Crops_Geotiff/", 
                      mcode[i], "_HarvAreaYield_Geotiff/", mcode[i], 
                      "_DataQuality_HarvestedArea.tif")
}
mdq <- rast(mfn_dq)

# function to select CV depending on Monf DQ
# heterogeneity = 0 for county level data
dqfun <- function(dq, cvc, cvs){
  ifelse(dq == 1, 0, 
         ifelse(dq %in% c(0.75, 0.5), cvs,
                ifelse(dq == 0.25, cvc, NA)))
}

# create DQI raster and save them to disk
outpath <- file.path(rd, "Monfreda/DataQualityIndex")
dir.create(outpath)

for(i in 1:length(mcode)){
  dqi <- lapp(c(mdq[[i]], rcvc, rcvs), fun = dqfun,
              filename = file.path(outpath, paste0(mcode[i], "_DQI.tif")), 
              overwrite = T, 
              wopt = list(names=mcode[i], filetype = "GTiff",
                          gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))
}

