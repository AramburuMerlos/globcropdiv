library(magrittr)
library(terra)
library(Recocrop)
library(data.table)

if(!grepl("globcropdiv$", getwd())){
  if (system('hostname', TRUE) %in% c("ESP-RH-9891", "LAPTOP-ST129J47")) { 
    setwd("G:/My Drive/globcropdiv/")
  } # else if { ... 
}

# Merge Data sets ####################################################
ecoc <- fread("AuxData/Ecocrops.csv")   # see Ecocrops_readme.txt 
monf <- fread("D:/Monfreda/Monfreda_Metadata.csv")
spam <- fread("D:/SPAM/SPAM_Metadata.csv")

## Ecocrop - Monfreda --------
d <- merge(monf, ecoc, by = "FAO_Code", all.x = T, all.y = F)
nas <- d[is.na(NAME),]
d <- d[!is.na(NAME),]
d2 <- merge(nas[,colnames(monf), with = F], ecoc, 
            by.x = "FAO_Code", by.y = "FAO_Code2", all.x = T, all.y = F)
nas <- d2[is.na(NAME),]
d2 <- d2[!is.na(NAME),]
setnames(d2, "FAO_Code.y", "FAO_Code2")
d <- rbind(d, d2)

d3 <- merge(nas[,colnames(monf), with = F], ecoc, 
            by.x = "FAO_Code", by.y = "FAO_Code3", all.x = T, all.y = F)
nas <- d3[is.na(NAME),] # see notes

fwrite(nas[,1:4], "AuxData/ExcMonfCrops.csv")

d3 <- d3[!is.na(NAME),]
setnames(d3, "FAO_Code.y", "FAO_Code3")
d <- rbind(d, d3)

# remove mushrooms because it is not a plant crop
d <- d[-which(Monf_Name == "mushroom"),]

rm(d2, d3)

## Ecocrop - SPAM -------
# change SPAM code for Coffee arabica and Panicum miliaceum  
d[SCIENTNAME == "Coffea arabica L.", SPAM_Code:= "acof"]
d[SCIENTNAME == "Panicum miliaceum L.", SPAM_Code:= "smil"]

# add SPAM names
d <- d[spam[,.SD, .SDcols = c("SPAM_Code", "SPAM_Name")], on = "SPAM_Code"]

# aggregated SPAM categories
todis <- c("ocer", "orts", "opul", "ooil", 
           "ofib", "trof", "temf", "vege", "rest")

# define data source in d. Add SPAM_Name. Order columns. Save it
d[, source:= ifelse(SPAM_Code %in% todis, "Monf", "SPAM")]
d[, crop:= ifelse(Monf_Name %in% c("coffee", "millet"), SPAM_Name, Monf_Name)]

# order columns and save
setcolorder(d, c("crop", "FAO_Name", "FAO_Code", "source", 
                 "SPAM_Name", "SPAM_Code", "Monf_Name"))
fwrite(d, "AuxData/CropSpecies.csv")

## Monfreda - SPAM -----
# Use Monfreda data for aggregated SPAM categories
# but following SPAM crop area, and avoiding Excluded Monfreda Crops. 

monf_crops <- d[source == "Monf", Monf_Name] %>% unique()
monf_sel <- monf[Monf_Name %in% monf_crops, ]
monf_sel[, source:= "Monf"]

spam_crops <- d[source == "SPAM", SPAM_Code] %>% unique()
spam_sel <- spam[SPAM_Code %in% spam_crops, ]
spam_sel[, source:= "SPAM"]
spam_sel[, Index:= NULL]

# add names
monf_sel <- merge(monf_sel, spam[,c("SPAM_Code", "SPAM_Name"), with = F], 
                  by = "SPAM_Code", all.x = T, all.y = F)

spam_sel <- merge(spam_sel, monf[,c("SPAM_Code", "Monf_Name"), with = F], 
                  by = "SPAM_Code", all.x = T, all.y = F)
spam_sel[SPAM_Code == "acof", Monf_Name:= "coffee"]
spam_sel[SPAM_Code == "smil", Monf_Name:= "millet"]

# columns order
cols.order <- c("Monf_Name", "source", "SPAM_Name", "SPAM_Code",  "FAO_Name",
                "FAO_Code",  "GROUP")
setcolorder(monf_sel, cols.order)
setcolorder(spam_sel, cols.order)

# combine and save
sel <- rbind(spam_sel, monf_sel)
sel[, crop:= ifelse(Monf_Name %in% c("coffee", "millet"), SPAM_Name, Monf_Name)]

setcolorder(sel, "crop")
fwrite(sel, "AuxData/CropCategories.csv")

# Total cropland ###############################################################
# Based on SPAM data
spam_abun <- "D:/SPAM/Physical_Area/*_A.tif" %>% Sys.glob() %>% rast()
totcl <- app(spam_abun, fun = sum, na.rm = T) 

# change 0 cropland to Na
totcl <- classify(totcl, rcl = cbind(0,NA))

# mask based on predictors to avoid future NA
clim <- "InData/WorldClim/2.1/wc5min/tavg/*.tif" %>% Sys.glob() %>% rast()
soil <- "InData/SoilGrids/phh2o/ph*_5min.tif" %>% Sys.glob() %>% rast()

totcl <- mask(totcl, clim[[1]])
totcl <- 
  mask(totcl, soil, 
       filename = "InData/TotalCropland.tif", overwrite = T, 
       wopt = list(names="TotalCropland", filetype = "GTiff",
                   gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))

# Mask SPAM Data ###############################################################
# based on total cropland
outdir <- "InData/CropAbundance"
dir.create(outdir, F, T)

sc <- sel[source == "SPAM", crop]

names(spam_abun) <- names(spam_abun) %>% 
  gsub("spam2010V2r0_global_A_", "", .) %>%
  gsub("_A", "", .) %>%
  tolower()

for(i in 1:length(sc)){
  scode <- sel[crop == sc[i], SPAM_Code]
  spam_abun[[scode]] <- 
    mask(spam_abun[[scode]], totcl,
         filename = file.path(outdir, paste0(sc[i], ".tif")), 
         overwrite = T, 
         wopt = list(names = sc[i], filetype = "GTiff",
                     gdal = c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6"))) 
}

# Fix Monfreda Data #######################################################
# readjust crop abundance data from Monfreda crops 
# so they add up to their corresponding SPAM crop. 
monf_abun <- "D:/Monfreda/Harv*175*/*/*_Har*AreaHectares.tif" %>% 
  Sys.glob() %>% rast()

names(monf_abun) <- names(monf_abun) %>%
  gsub("_HarvestedAreaHectares", "", .) 

for(i in 1:length(todis)){
  # SPAM values
  v <- values(spam_abun[todis[i]], dataframe = TRUE) %>% setDT()
  setnames(v, names(v), "stot")
  v[, cell:= 1:.N]
  v <- v[!is.na(stot)]
  
  # Monfreda crops within this spam category
  mc <- sel[SPAM_Code == todis[i], Monf_Name]
  v[, (mc):= extract(monf_abun[mc], cell)]
  
  # overall proportion of Monfreda crops (for cells with no Monf data)
  prop <- colSums(v[,..mc], na.rm = T) %>% `/`(sum(.))
  
  # Monf totals
  v[, mtot:= rowSums(.SD, na.rm = T), .SDcols = mc]
  
  # Update Monf Values:
  ## 1. Change all Monf crops to 0 when stot is 0
  v[stot == 0, (mc):= 0]
  
  ## 2. Recalculate Monf crops abundance to match spam tot
  v[mtot>0 & stot>0, (mc):= lapply(.SD, function(x) x*stot/mtot), .SDcols = mc]
  
  ## 3. Use total crop proportinos for cell where stot >0 but no monf data
  v[mtot == 0 & stot > 0, (mc):= lapply(prop, function(x) x * stot)]

  # set values and write raster
  for(j in 1:length(mc)){
    r <- rast(monf_abun[[mc[j]]])
    vj <- rep(NA_real_, ncell(r))
    vj[v$cell] <- v[[mc[j]]]
    values(r) <- vj
    wopt <- list(names = mc[j], filetype = "GTiff",
                 gdal = c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6"))
    writeRaster(r, filename = file.path(outdir, paste0(mc[j], ".tif")), 
                overwrite = T, wopt = wopt) 
  }
}

# Notes ---------

# Crop categories in Monfreda data that aren't represented in Ecocrops
# i. Species not included in Ecocrops data
# - Triticale
# ii. Categories with no specific/explicit species
# - Mixed grain
# - Oil Seeds for forage
# - Forage Products nes
# - Vegetables, roots fodder nes

# Note that some categories may include species not available in Ecocrops.
# (particularly "nes" categories) 
# Ecocrop species may be repeated if they are included in more than one Monfreda category
# Monfreda categories may be repeated if they include more than one crop species





