library(magrittr)
library(data.table)
library(terra)

if(!getwd() %like%  "globcropdiv$") warning("See 0000_wd.R")

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

