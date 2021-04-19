library(magrittr)
library(data.table)
library(terra)

if(!grepl("globcropdiv$", getwd())){
  if (system('hostname', TRUE) %in% c("ESP-RH-9891", "LAPTOP-ST129J47")) { 
    setwd("G:/My Drive/globcropdiv/")
  } # else if { ... 
}


# Maxent suitability
mxt <- "OutData/SDM/*.tif" %>%  Sys.glob() %>% rast() 

# Ecocrop Suitability 
eco <- "OutData/Ecocrop/*.tif" %>% Sys.glob() %>% rast()

outpath <- "OutData/Int_Suit"
dir.create(outpath, F, T)

all.equal(names(mxt), names(eco))
crops <- names(eco)

for(i in 1:length(crops)){
  app(c(eco[[i]], mxt[[i]]), fun = "min",
       filename = file.path(outpath, paste0(crops[i], ".tif")), 
       overwrite = T, 
       wopt = list(names = crops[i], filetype = "GTiff",
                   gdal = c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))
}

