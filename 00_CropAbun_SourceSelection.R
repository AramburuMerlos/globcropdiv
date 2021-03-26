# install.packages(c("meteor", "terra", "here", "data.table"))
# remotes::install_github("cropmodels/Recocrop")
library(magrittr)
library(terra)
library(Recocrop)
library(data.table)

dirn <- "AuxData"
dir.create(dirn, F)

d <- "parameters/ecocrop.rds" %>% system.file(package="Recocrop") %>% readRDS
setDT(d)
d[, NAME:= sub(" \\*", "", NAME)]
d[, NAME:= sub("\\*", "", NAME)]
fwrite(d, file.path(dirn, "EcoCropPars.csv"))
rm(d)
# FAO_Codes and categories added externally
# see Ecocrops_metadata.txt for a description of the following table

# Merge Data sets ####################################################
# Ecocrop - Monfreda - SPAM
# Define selected crop categories and their corresponding Ecocrop species

# Match EcoCrops with Monfreda and SPAM data (~FAO) crop groups
ecoc <- fread(file.path(dirn, "Ecocrops.csv"))

monf <- fread("D:/Monfreda/Cropcat.csv")
monf <- monf[!is.na(FAO_Code),]

spam <- fread("D:/SPAM/4-Methodology-Crops-of-SPAM-2005-2015-02-26.csv")

spam[,FAO_Code:= as.numeric(FAOCODE)] #NA for categ with >1 FAO cat
spam <- spam[!is.na(FAO_Code),]
setnames(spam, 
         c("SPAM short name", "SPAM long name", "FAONAMES"),
         c("SPAM_Code","SPAM_Name", "FAO_Name"))
spam[, c("No. crt. ", "FAOCODE"):= NULL]

# does all SPAM categories match at least one Monf category
all(spam$FAO_Code %in% monf$FAO_Code)

# is there a SPAM category that matches more than one Monf Category?
max(sapply(spam$FAO_Code, function(x)sum(x %in% monf$FAO_Code)))

# remove Monfreda categories that are available in SPAM
mns <- monf[!monf$FAO_Code %in% spam$FAO_Code,]

# combine (bind) spam and monfreda crop categories
spmo <- rbind(mns, spam, fill = T)

# merge cobmined crop categories with ecocrop database
d <- merge(spmo, ecoc, by = "FAO_Code", all.x = T, all.y = F)
nas <- d[is.na(NAME),]
d <- d[!is.na(NAME),]
d2 <- merge(nas[,colnames(spmo), with = F], ecoc, 
            by.x = "FAO_Code", by.y = "FAO_Code2", all.x = T, all.y = F)
nas <- d2[is.na(NAME),]
d2 <- d2[!is.na(NAME),]
setnames(d2, "FAO_Code.y", "FAO_Code2")
d <- rbind(d, d2)

d3 <- merge(nas[,colnames(spmo), with = F], ecoc, 
            by.x = "FAO_Code", by.y = "FAO_Code3", all.x = T, all.y = F)
nas <- d3[is.na(NAME),] # see notes

fwrite(nas[,1:4], file.path(dirn, "ExcMonfCrops.csv"))


d3 <- d3[!is.na(NAME),]
setnames(d3, "FAO_Code.y", "FAO_Code3")
d <- rbind(d, d3)

d[,source:= ifelse(is.na(SPAM_Code), "Monf", "SPAM")]
setcolorder(d, c("FAO_Name", "FAO_Code", "source", 
                 "SPAM_Name", "SPAM_Code", "Monf_Name"))

# remove mushrooms because it is not a plant crop
d <- d[-which(Monf_Name == "mushroom"),]

fwrite(d, file.path(dirn, "EcocropSpeciesWithAbundanceData.csv"))

# Selected crops ############

spmo[,source:= ifelse(is.na(SPAM_Code), "Monf", "SPAM")]
spmo[, crop:= ifelse(is.na(SPAM_Code), Monf_Name, SPAM_Name)]
setcolorder(spmo, c("crop", "FAO_Name", "FAO_Code", "source",  
                    "SPAM_Name", "SPAM_Code", "Monf_Name"))
# remove mushrooms because it is not a plant crop
spmo <- spmo[-which(Monf_Name == "mushroom"),]

# add file.path and names
spmo[, file := ifelse(source == "Monf",
                      paste0("Monfreda/HarvestedAreaYield175Crops_Geotiff/", 
                             Monf_Name, "_HarvAreaYield_Geotiff/", 
                             Monf_Name, "_HarvestedAreaHectares.tif"),
                      paste0("SPAM/Physical_Area/spam2010V2r0_global_A_",
                             toupper(SPAM_Code),"_A.tif"))]

fwrite(spmo, file.path(dirn, "CropAbundanceSource.csv"))

# Total Cropland (ha per cell) --------
allcrops <- rast(file.path("D:", spmo$file))
totcl <- app(allcrops, fun = sum, na.rm = T) 
classify(totcl, rcl = cbind(0,NA), 
         filename = "InData/TotalCropland.tiff", overwrite = T, 
         wopt = list(names="TotalCropland", filetype = "GTiff",
                     gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))

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





