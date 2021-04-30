library(magrittr)
library(data.table)
library(terra)

if(!grepl("globcropdiv$", getwd())){
  if (system('hostname', TRUE) %in% c("ESP-RH-9891", "LAPTOP-ST129J47")) { 
    setwd("G:/My Drive/globcropdiv/")
  } # else if { ... 
}

source("Functions/allocate.R")

# Data prep -------

d <- "InData/TotalCropland.tif" %>%
  rast() %>%
  values() %>% 
  is.na() %>% `!` %>% 
  which() %>% 
  data.table(cell = .)

# Suitability
r_suit <- "OutData/SDM/*_AVG.tif" %>% Sys.glob() %>% rast()
crops <- names(r_suit)
scols <- paste0("s.", crops)
d[, (scols):= extract(r_suit, cell)]
 
# total cropland per cell
d[, tcl:= extract(rast("InData/TotalCropland.tif"), cell)]

# total crop area per crop
r_area <- "InData/CropAbundance/*.tif" %>% Sys.glob() %>% rast()
acols <- paste0("a.", crops)
d[, (acols):= extract(r_area, cell)]
crop_area <- colSums(d[, ..acols])
names(crop_area) <- crops
d[, (acols):= NULL]


#################
d <- allocate(d, crops, crop_area,  
              scols = paste0("s.", crops),
              tccol = "tcl")

fwrite(d, "OutData/allocated_sdm.csv")


