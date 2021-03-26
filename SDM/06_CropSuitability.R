library(data.table)
library(terra)
library(maxnet)

if(!getwd() %like%  "globcropdiv$") warning("See 0000_wd.R")

# directory of input data
indir <- "D:"

# functions ----
source("Functions/maxent_output_transformation.R")
source("Functions/sample_presences.R")

# DATA PREP ###########################################################
## Predictors --------------
preds.fn <- fread("AuxData/SelectedPreds.csv")$file_name
preds <- rast(file.path(indir, preds.fn))
names(preds) <- fread("AuxData/SelectedPreds.csv")$short_name

## Total Cropland (per cell) --------
# predictions will be made for cropland only
totcl <- rast("InData/TotalCropland.tiff")

# Mask Total Cropland based on predictors
pred_mask <- app(preds, sum, na.rm = FALSE) # to get NA if any layer has NA
totcl <- mask(totcl, pred_mask)
rm(pred_mask)

## create data.table ------------
d <- data.table(cell = 1:ncell(totcl), tot = values(totcl)[,1])

# remove NA (non cropland) cells, order by cell number
d <- d[!is.na(tot),]

# Extract crop abundances
sc <- fread("AuxData/CropAbundanceSource.csv")
rabun <- rast(file.path(indir, sc$file))
crops <- sc$crop
names(rabun) <- crops

d[, (crops):= extract(rabun, cell)]

# set NA abundances to 0
for(j in crops) set(d, i = d[,.I[is.na(.SD)], .SDcols = j], j = j, 0)

# change abundances to fraction
d[, (crops):= lapply(.SD, function(x) x/tot), .SDcols = (crops)]

## Extract predictor values --------
preds_names <- names(preds)
d[, (preds_names):= extract(preds, cell)]


# output dir 
outdir <- "OutData/SDM/Suitability"
dir.create(outdir, F, T)


# RUN MODELS  ############################################################

# maxent parameters from regional CV
j_pm = 2 ; regmult = 3.5 ; ss = 0.04 ; npr = 4000

tr <- which(d$cell%%(1/ss) == 0)

set.seed(1)

for(i in 1:length(crops)){
  pr <- fpr(d[tr], crops[i], j_pm, npr)
  m <- maxnet(c(rep(1, npr), rep(0, d[tr][,.N])),
              d[tr][c(pr, 1:length(tr)), ..preds_names],
              regmult = regmult)  
  logp <- predict(m, d[, ..preds_names], type = "logistic")
  d[, paste0(crops[i],".suit"):= logp]
  rm(pr, m, logp); gc()
  print(paste(i, "/", length(crops)))
}

### Write on Disk ----
for(i in 1:length(crops)){
  r <- rast(totcl)
  v <- rep(NA_real_, ncell(r))
  isuit <- paste0(crops[i], ".suit")
  v[d$cell] <- d[[isuit]]
  values(r) <- v
  fn <- file.path(outdir, paste0(isuit, ".tif"))
  wopts <- list(names = isuit, filetype = "GTiff",
                gdal=c("COMPRESS=Deflate","PREDICTOR=1", "ZLEVEL=6"))
  writeRaster(r, filename = fn, overwrite = T, wopt = wopts)
}
