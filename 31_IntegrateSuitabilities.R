library(magrittr)
library(terra)

if(system('hostname', TRUE) == "ESP-RH-9891"){
  setwd("D:/globcropdiv/")
} else if(system('hostname', TRUE) == "LAPTOP-ST129J47"){ 
  setwd("G:/My Drive/globcropdiv/")
} # else if { ... 

# SDM model averaged suitability
sdm <- "OutData/SDM/*_AVG.tif" %>%  Sys.glob() %>% rast() 

# Ecocrop Suitability 
eco <- "OutData/Ecocrop/*.tif" %>% Sys.glob() %>% rast()

outpath <- "OutData/Int_Suit"
dir.create(outpath, F, T)

all.equal(names(sdm), names(eco))
crops <- names(eco)

for(i in 1:length(crops)){
  app(c(eco[[i]], sdm[[i]]), fun = "min",
       filename = file.path(outpath, paste0(crops[i], ".tif")), 
       overwrite = T, 
       wopt = list(names = crops[i], filetype = "GTiff",
                   gdal = c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))
}

