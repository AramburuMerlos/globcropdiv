library(magrittr)
library(data.table)
library(terra)

if(!grepl("globcropdiv$", getwd())){
  if (system('hostname', TRUE) %in% c("ESP-RH-9891", "LAPTOP-ST129J47")) { 
    setwd("G:/My Drive/globcropdiv/")
  } # else if { ... 
}


# Maxent suitability
mxt <- "OutData/SDM/Suitability/*.tif" %>%  Sys.glob %>% rast 


# Ecocrop Suitability 
eco <- "OutData/EcocropSuit/byCrop/*.tif" %>% Sys.glob %>% rast


outpath <- "OutData/Int_Suit"
dir.create(outpath, F, T)


# match names
imxt <- names(mxt) %>% gsub(".suit", "", .) %>% match(names(eco), .)


for(i in 1:length(imxt)){
  app(c(eco[[i]], mxt[[imxt[i]]]), fun = "min",
       filename = file.path(outpath, paste0(names(eco)[i], ".tif")), 
       overwrite = T, 
       wopt = list(names = names(eco)[i], filetype = "GTiff",
                   gdal = c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))
}

