library(data.table)
library(raster)
library(terra)
library(rnaturalearth)
library(rnaturalearthhires)

nwd <- nchar(getwd())
if(substr(getwd(),  nwd - 10, nwd) != "globcropdiv") stop("See 0000_wd.R")

sc <- fread("AuxData/CropAbundanceSource.csv")

# Drive or folder where the Raw SPAM and Monfreda data is saved
rd <- "D:/"

# Filter Monfreda Data ########################################
# filter out crop abundance cells that come from "very aggregated data" (VAD)
# VAD is defined as a function of the CV for annual rainfall (Bio12) and GDD
# within the cropland in which the crop abundance data is spread
# Procedure: get max(CV T, CVrain) within each polygon for each aggregation level, 
# considering only the cropland of each polygon
# then filter data that comes from polygons with more than X CV (30%)


# upload Monfreda crop abundance data
mcode <- sc[source == "Monf", Monf_Name]

mfn <- vector(length = length(mcode))
for(i in 1:length(mcode)){
  mfn[i] <- paste0(rd, "Monfreda/HarvestedAreaYield175Crops_Geotiff/", 
               mcode[i], "_HarvAreaYield_Geotiff/", mcode[i], 
               "_HarvestedAreaHectares.tif")
}
mcrops <- rast(mfn)
  
# rast with total cropland (ha) per cell based on Monfreda data
monf_tot <- app(mcrops, fun = sum, na.rm = T)

# raster with 1 for cropland and NA for non cropland
monf_cropland <- monf_tot
values(monf_cropland) <- values(monf_tot, mat = FALSE) > 0

# Country CV Mask ----
# rast with countries 
spcountries <- ne_countries()
spcountries$ID <- 1:nrow(spcountries)
rcountries <- rasterize(vect(spcountries), mcrops, field = "ID") 

# Growing Degree Days and Precipitation
rain <- rast(file.path(rd, "WorldClim/2.1/wc5min/bioc/wc2.1_5m_bio_12.tif"))
GDD <- rast(file.path(rd, "WorldClim/2.1/wc5min/extra/GDD.tif"))

# mask cropland
rain <- mask(rain, monf_cropland, maskvalue = 0)
GDD <- mask(GDD, monf_cropland, maskvalue = 0)
# for rain, mask irrigated land
irr <- rast(paste0(rd, "AQUASTAT/gmia_v5_aei_pct.asc"))
irr_m <- irr > 50
rain <- mask(rain, irr_m, maskvalue = 1)
names(rain) <- "rain"

gdd_country <- crosstab(raster(rcountries), 
                        raster(GDD), 
                        long = TRUE, digits = 0)

rain_country <- crosstab(raster(rcountries), 
                         raster(rain), 
                         long = TRUE, digits = 0)

# add 100mm to all annual precipitation data to avoid high CV because of low Pp
rain_country$rain <- rain_country$rain + 100

# function to get Cv from freq table
cvfun <- function(x, f){
  tf <- sum(f)
  avg <-  sum(x * f/tf)
  var <- sum((x - avg)^2 * f/tf)
  return(sqrt(var)/avg)
}

setDT(rain_country)
CVrain_country <- rain_country[, .(CVrain = cvfun(rain, Freq)), by = .(ID)]

setDT(gdd_country)
CVgdd_country <- gdd_country[, .(CVgdd = cvfun(GDD, Freq)), by = .(ID)]

CV_country <- merge.data.table(CVrain_country, CVgdd_country, by = "ID")
CV_country[, maxCV := pmax(CVrain, CVgdd)]

rcm_country <- as.matrix(CV_country[,.(ID,maxCV)])

# raster with max CV (between rainfall and temperature) per country (across their cropland)
rcvc <- classify(rcountries, rcm_country, othersNA = TRUE)
country_mask <- rcvc
values(country_mask) <- values(rcvc, mat = FALSE) > 0.30


# States CV Mask ----
# rast with states / provinces
spstates <- ne_states()
spstates$ID <- 1:nrow(spstates)
rstates <- rasterize(vect(spstates), mcrops, field = "ID") 
gdd_states <- crosstab(raster(rstates), 
                       raster(GDD), 
                       long = TRUE, digits = 0)

rain_states <- crosstab(raster(rstates), 
                        raster(rain), 
                        long = TRUE, digits = 0)

# add 100mm to all annual precipitation data to avoid high CV because of low Pp
rain_states$rain <- rain_states$rain + 100

setDT(rain_states)
CVrain_states <- rain_states[, .(CVrain = cvfun(rain, Freq)), by = .(ID)]

setDT(gdd_states)
CVgdd_states <- gdd_states[, .(CVgdd = cvfun(GDD, Freq)), by = .(ID)]

CV_states <- merge.data.table(CVrain_states, CVgdd_states, by = "ID")
CV_states[, maxCV := pmax(CVrain, CVgdd)]

rcm_states <- as.matrix(CV_states[,.(ID,maxCV)])

# raster with max CV (between rainfall and temperature) per country (across their cropland)
rcvs <- classify(rstates, rcm_states, othersNA = TRUE)
state_mask <- rcvs
values(state_mask) <- values(rcvs, mat = FALSE) > 0.3

# Monfreda data quality
mfn_dq <- vector(length = length(mcode))
for(i in 1:length(mcode)){
  mfn_dq[i] <- paste0(rd, "Monfreda/HarvestedAreaYield175Crops_Geotiff/", 
                      mcode[i], "_HarvAreaYield_Geotiff/", mcode[i], 
                      "_DataQuality_HarvestedArea.tif")
}
dq_mcrops <- rast(mfn_dq)

# mask data from aggregated and heterogeneous sources and save it to disk
out.dir <- file.path(rd, "WorldCropAbundance")
dir.create(out.dir)

out.dir.ha <- file.path(rd, "WorldCropAbundance/hectares")
dir.create(out.dir.ha)

for(i in 1:length(mcode)){
  # change values lower than 0.01 ha - 100 m2; 10m by 10m - to NA
  mc <- classify(mcrops[[i]], cbind(0,0.01,NA), include.lowest = TRUE)
  # create crop specific heterogeneous country/state data mask
  cm <- sm <- dq_mcrops[[i]]
  vdq <- values(dq_mcrops[[i]], mat=FALSE)
  values(cm) <- vdq == 0.25 & values(country_mask, mat=FALSE) == 1
  values(sm) <- (vdq == 0.5 | vdq == 0.75) & values(state_mask, mat=FALSE) == 1
  im <- max(cm, sm)
  outfn <- file.path(out.dir.ha, paste0(mcode[i], "_MONF_ha.tif"))
  mask(mc, im, maskvalue = 1,
       filename = outfn, 
       overwrite = T, 
       wopt = list(names = mcode[i], filetype = "GTiff",
                   gdal = c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))
}


# Filter SPAM Data ###########################################

# upload SPAM data
scode <- sc[source == "SPAM", SPAM_Code]
sfn <- vector(length = length(scode))
for(i in 1:length(scode)){
  sfn[i] <- paste0(rd, "SPAM/Physical_Area/spam2010V2r0_global_A_", 
                   toupper(scode[i]), "_A.tif")
}
scrops <- rast(sfn)
compareGeom(scrops, mcrops, lyrs=FALSE)

# remove values lower than 0.01 ha and save to disk
for(i in 1:length(scode)){
  outfn <- file.path(out.dir.ha, paste0(scode[i], "_SPAM_ha.tif"))
  classify(scrops[[i]], cbind(0,0.01,NA), include.lowest = TRUE,
           filename = outfn, 
           overwrite = T, 
           wopt = list(names = scode[i], filetype = "GTiff",
                       gdal = c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6"))) 
}

crops <- c(scrops, mcrops)

# create cropland mask with total cropland (in ha)
total0 <- app(crops, fun = sum, na.rm = T)
total <- classify(total0, cbind(0,NA),
                  filename = file.path(out.dir, "Total_Cropland_ha.tif"),
                  overwrite = T, 
                  wopt = list(names = "Total Cropland", filetype = "GTiff"))


# Area to Cropland Fraction #################
fn <- Sys.glob(file.path(out.dir.ha,"*.tif"))
area_ha <- rast(fn)
codes <- gsub(paste0(out.dir.ha,"/"), "", gsub("_ha.tif$", "", fn))
total <- rast(file.path(out.dir, "Total_Cropland_ha.tif"))

out.dir.fr <- file.path(rd, "WorldCropAbundance/fraction")
dir.create(out.dir.fr)
vtotal <- values(total, mat = FALSE)

for(i in 1:nlyr(area_ha)){
  outfn <- file.path(out.dir.fr, paste0(codes[i], "_fr.tif"))
  outr <- area_ha[[i]]
  values(outr) <- values(area_ha[[i]], mat = FALSE)/vtotal
  writeRaster(outr,
              filename = outfn, 
              overwrite = T, 
              wopt = list(names = codes[i], filetype = "GTiff",
                          gdal = c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6"))) 
}

