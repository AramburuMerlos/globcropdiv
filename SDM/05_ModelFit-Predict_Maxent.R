# the probability of deining a given cell as presence will be in function of its crop fraction. (p(cell = 1) = crop_fraction)
# then, the resulting absence cells will be subsampled to reduce the amount of training data
# here, the sampleing effort is the number of trials in each cell (size in rbinom)

library(data.table)
library(terra)
library(maxnet)

nwd <- nchar(getwd())
if(substr(getwd(),  nwd - 10, nwd) != "globcropdiv") warning("See 0000_wd.R")

# directory of input data
indir <- "D:"
cpath <- file.path(indir, "WorldCropAbundance")

# Crop Abundance Data ----------
sc <- fread("AuxData/CropAbundanceSource_Filtered.csv")
# KEEP ONLY CROPS WITH >75% Of THEIR DATA NON-FILTERED 
sel_sc <- sc[nf_prop_ha > 0.75,]
crops <- sel_sc[, filename]
crops.fn <- paste0(crops, "_", sel_sc[, source])
fn <- file.path(cpath, paste0("fraction/", crops.fn, "_fr.tif"))
abund <- rast(fn)

# Predictors --------------
preds.fn <- fread("AuxData/SelectedPreds.csv")$file_name
preds <- rast(file.path(indir, preds.fn))
names(preds) <- fread("AuxData/SelectedPreds.csv")$short_name

# compare predictors and abundance Geometries
compareGeom(preds, abund)

# Cropland Mask --------
cmask <- rast(file.path(cpath, "Total_Cropland_ha.tif"))
preds <- mask(preds, cmask)

# Mask cropland-mask based on predictors
pred_mask <- app(preds, sum, na.rm = FALSE) # to get NA if any layer has NA
cmask <- mask(cmask, pred_mask)

# Predictors data.table to fit model
dpreds <- data.table(cell = seq_len(ncell(preds)), values(preds))
dpreds <- dpreds[complete.cases(dpreds),]

# number of absence cells: 
nab <- 2e4

# Set sampling efforts ------
# sampling effort is number of trials in each cell
# ntrials is a function of total crop area. Rare crops -> more trials
arealogs <- round(log(sel_sc[, nonFiltered_ha]))
ntrials <- max(arealogs) - arealogs + 1

# run MAXENT -----
outpath <- file.path("OutData/SDM/Maxent")
dir.create(outpath)

r <- rast(cmask)
# loop for all crops  ----
set.seed(1234)
np <- vector(length = length(crops)) # vector to save number of presence values
for(i in 1:length(crops)){
  d <- data.table(cell = (1:ncell(abund[[i]])),
                  frac = values(abund[[i]])[,1])
  # (un)mask cells using crop mask, add 0 in cropland cells without crop (NA in abund[[i]])
  d[,cm:= values(cmask)]
  d <- d[!is.na(cm),]
  d[,cm:= NULL]
  d[, frac:= ifelse(is.na(frac), 0, frac)]
  d[, frac:= pmin(frac, 1)] # limit fraction to a max of 1 
  
  # is the crop detected in each cell?
  d[, pres := rbinom(.N, ntrials[i], frac)]
  # presence cells 
  pcells <- d[pres > 0, cell]
  # absence cells (keep only nab absence observations)
  acells <- sample(d[pres == 0, cell], nab)
  
  # keep track of number of presence cells 
  np[i] <- length(pcells)
  
  # predictor for absences
  ap <- as.data.table(extract(preds, acells))
  # predictors for presence cells
  pp <- as.data.table(extract(preds, pcells))
  # create data table with predictors and vector with 0/1 values
  d <- rbindlist(list(ap, pp))
  p <- c(rep(0L, nrow(ap)), rep(1L, nrow(pp)))
  
  # fit
  m <- maxnet(p, d)
  
  # Predict, set values in raster, save
  mp <- predict(m, dpreds[,-1], type = "cloglog") 
  v <- rep(NA_real_, ncell(preds))
  v[dpreds$cell] <- mp[,1]
  values(r) <- v
  writeRaster(r,   
              filename = file.path(outpath, paste0(crops[i], "_Maxent_pred.tif")),
              overwrite = T, 
              wopt = list(names=crops[i], filetype = "GTiff",
                          gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))
}
  
sel_sc[, npres_Maxent:= np]
sel_sc[, ntrials_Maxent:= ntrials]
fwrite(sel_sc, "AuxData/CroplandSampling_Maxent.csv")
plot(log(npres_Maxent) ~ log(nonFiltered_ha), data = sel_sc, 
     col = rainbow(max(ntrials))[ntrials], pch = 19)
