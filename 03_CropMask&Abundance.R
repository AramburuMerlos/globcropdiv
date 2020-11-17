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

# Directories
out.dir <- file.path(rd, "WorldCropAbundance")
dir.create(out.dir)

out.dir.ha <- file.path(rd, "WorldCropAbundance/hectares")
dir.create(out.dir.ha)

out.dir.fr <- file.path(rd, "WorldCropAbundance/fraction")
dir.create(out.dir.fr)


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
monf_cropland <- arith(monf_tot, fun = function(x) x > 0)

# Country and State CV Mask ----
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
irr <- rast(paste0(rd, "AQUASTAT/gmia_v5_aei_pct.asc"))
irr_m <- arith(irr, function(x) x > 50)
rain <- mask(rain, irr_m, maskvalue = 1)
names(rain) <- "rain"

DT <- data.table(values(c(rcountries, rstates, GDD, rain)))
DT <- DT[complete.cases(DT),]

# add 100mm to all annual precipitation data to avoid high CV because of low Pp
DT[, rain:= rain + 100]

# compute CV by country, get max CV and create mask
cvf <- function(x)sd(x, na.rm = T)/mean(x, na.rm = T)
CVc <- DT[, lapply(.SD, cvf), .SDcols = c("GDD", "rain"), by = "cID"]
CVc[, maxCV := pmax(GDD, rain)]
rcm_country <- as.matrix(CVc[,.(cID,maxCV)])
rcvc <- classify(rcountries, rcm_country, othersNA = TRUE)
country_mask <- arith(rcvc, function(x) x > 0.30)

# compute CV by state, get max CV and create mask
CVs <- DT[, lapply(.SD, cvf), .SDcols = c("GDD", "rain"), by = "sID"]
CVs[, maxCV := pmax(GDD, rain)]
rcm_state <- as.matrix(CVs[,.(sID, maxCV)])
rcvs <- classify(rstates, rcm_state, othersNA = TRUE)
state_mask <- arith(rcvs, function(x) x > 0.30)

# Monfreda data quality
mfn_dq <- vector(length = length(mcode))
for(i in 1:length(mcode)){
  mfn_dq[i] <- paste0(rd, "Monfreda/HarvestedAreaYield175Crops_Geotiff/", 
                      mcode[i], "_HarvAreaYield_Geotiff/", mcode[i], 
                      "_DataQuality_HarvestedArea.tif")
}
dq_mcrops <- rast(mfn_dq)

# create crop specific heterogeneous country/state data mask
maskfun <- function(dq, cm, sm){
  (dq == 0.25 & cm == 1) | ((dq == 0.5 | dq == 0.75) & sm == 1)
}

# mask data from aggregated and heterogeneous sources and save it to disk

for(i in 1:length(mcode)){
  # change values lower than 0.01 ha (i.e. 100 m2; 10m by 10m) to NA
  mc <- classify(mcrops[[i]], cbind(0,0.01,NA), include.lowest = TRUE)
  # create crop specific mask and save masked crop abundance to folder
  mi <- lapp(c(dq_mcrops[[i]], country_mask, state_mask), fun = maskfun)
  outfn <- file.path(out.dir.ha, paste0(mcode[i], "_MONF_ha.tif"))
  outr <- mask(mc, mi, maskvalue = 1,
               filename = outfn, overwrite = T, 
               wopt = list(names=mcode[i], filetype = "GTiff",
                           gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))
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

# Total Cropland #########################

## compute total number of ha and cells with data for original and masked crop abundance
# Total hectares Monf crops
thm <- global(mcrops, sum, na.rm = T) 
thm <- data.table(fullname = rownames(thm), Total_ha = thm$sum)
thm[, Monf_Name:= gsub("_HarvestedAreaHectares", "", fullname)]
thm[, fullname:= NULL]
# Total # "presence" cells Monf crops
tpm <- global(mcrops, fun = function(x)sum(x > 0, na.rm = T)) 
tpm <- data.table(fullname = rownames(tpm), Total_p = tpm$global)
tpm[, Monf_Name:= gsub("_HarvestedAreaHectares", "", fullname)]
tpm[, fullname:= NULL]

# Total hectares SPAM crops
ths <- global(scrops, sum, na.rm = T) 
ths <- data.table(fullcode = rownames(ths), Total_ha = ths$sum)
ths[, SPAM_Code:= tolower(substr(fullcode, 23,26))]
ths[, fullcode:= NULL]

# Total # "presence" cells SPAM crops
tps <- global(scrops, fun = function(x)sum(x > 0, na.rm = T)) 
tps <- data.table(fullcode = rownames(tps), Total_p = tps$global)
tps[, SPAM_Code:= tolower(substr(fullcode, 23,26))]
tps[, fullcode:= NULL]

# Total hectares and # of presence in filtered/masked data
fh <- global(crops, fun = function(x) sum(x, na.rm = T)) 
fp <- global(crops, fun = function(x) sum(x > 0, na.rm = T)) 
df <- data.table(name = names(crops), 
                 nonFiltered_ha = fh$global, 
                 nonFiltered_p = fp$global)

# merge results with source summary data. 
sc <- merge.data.table(sc, thm, by = "Monf_Name", all.x = T)
sc <- merge.data.table(sc, tpm, by = "Monf_Name", all.x = T)
sc <- merge.data.table(sc, ths, by = "SPAM_Code", all.x = T)
sc <- merge.data.table(sc, tps, by = "SPAM_Code", all.x = T)

# collapse columns
sc[,Total_ha:= pmax(Total_ha.x, Total_ha.y, na.rm = T)]
sc[,Total_p:= pmax(Total_p.x, Total_p.y, na.rm = T)]
sc[, `:=`(Total_p.x = NULL, Total_p.y = NULL, Total_ha.x = NULL, Total_ha.y = NULL)]

## Hectares and number of presence for Filtered data 
# upload all filtered crop data and compute totals
crops <- rast(Sys.glob(file.path(out.dir.ha, "*.tif")))
thc <- global(crops, function(x) sum(x, na.rm = T))
thc <- data.table(filename = rownames(thc), nonFiltered_ha = thc$global)

tpc <- global(crops, function(x) sum(x > 0, na.rm = T))
tpc <- data.table(filename = rownames(tpc), nonFiltered_p = tpc$global)

# merge results
sc[, filename:= ifelse(nchar(SPAM_Code) == 0, Monf_Name, SPAM_Code)]
sc <- merge.data.table(sc, thc, by = "filename", all.x = T)
sc <- merge.data.table(sc, tpc, by = "filename", all.x = T)
sc[, nf_prop_ha:= nonFiltered_ha/Total_ha]
sc[, nf_prop_p:= nonFiltered_p/Total_p]

fwrite(sc, "AuxData/CropAbundanceSource.csv")


# create cropland mask with total cropland (in ha)
crops <- rast(Sys.glob(file.path(out.dir.ha, "*.tif")))
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

for(i in 1:nlyr(area_ha)){
  outfn <- file.path(out.dir.fr, paste0(codes[i], "_fr.tif"))
  lapp(c(area_ha[[i]], total), function(x, y) x/y,
       filename = outfn, 
       overwrite = T, 
       wopt = list(names = codes[i], filetype = "GTiff",
                   gdal = c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6"))) 
}

