# install.packages(c("meteor", "terra", "here", "data.table"))
# remotes::install_github("cropmodels/Recocrop")
library(terra)
library(Recocrop)
library(data.table)

dirn <- "AuxData"
dir.create(dirn)

fname <- system.file("parameters/ecocrop.rds", package="Recocrop")
d <- readRDS(fname)
setDT(d)
d[, NAME:= sub(" \\*", "", NAME)]
d[, NAME:= sub("\\*", "", NAME)]
fwrite(d, file.path(dirn, "EcoCropPars.csv"))
rm(fname, d)
# FAO_Codes and categories added externally
# see Ecocrops_metadata.txt for a description of the following table

# Merge Data sets ####################################################

# Match EcoCrops with Monfreda data (~FAO) crop groups
ecrops <- fread(file.path(dirn, "Ecocrops.csv"))

if(dir.exists("D:/Monfreda")) {
  monfcrops <- fread("D:/Monfreda/CropCat.csv")  
} else {
  monfcrops <- fread("InData/Monfreda/CropCat.csv")
}

setnames(monfcrops, "CROPNAME", "Monf_Name")
d <- merge(monfcrops, ecrops, by = "FAO_Code", all.x = T, all.y = F)
nas <- d[is.na(NAME),]
d <- d[!is.na(NAME),]
d2 <- merge(nas[,colnames(monfcrops), with = F], ecrops, 
            by.x = "FAO_Code", by.y = "FAO_Code2", all.x = T, all.y = F)
nas <- d2[is.na(NAME),]
d2 <- d2[!is.na(NAME),]
setnames(d2, "FAO_Code.y", "FAO_Code2")
d <- rbind(d, d2)

d3 <- merge(nas[,colnames(monfcrops), with = F], ecrops, 
            by.x = "FAO_Code", by.y = "FAO_Code3", all.x = T, all.y = F)
nas <- d3[is.na(NAME),] # see notes
fwrite(nas[,1:4], file.path(dirn, "ExcMonfCrops.csv"))


d3 <- d3[!is.na(NAME),]
setnames(d3, "FAO_Code.y", "FAO_Code3")
d <- rbind(d, d3)
fwrite(d, file.path(dirn, "Monfcrops.csv"))


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





