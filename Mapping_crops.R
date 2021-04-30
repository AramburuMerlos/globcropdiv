library(magrittr)
library(data.table)
library(terra)
library(rnaturalearth)
library(colorRamps)

if(!grepl("globcropdiv$", getwd())){
  if(system('hostname', TRUE) %in% c("ESP-RH-9891", "LAPTOP-ST129J47")){ 
    setwd("G:/My Drive/globcropdiv/")
  } # else if { ... 
}

countries <- ne_download(scale = 10, type = "countries")

# Actual crop area --------
# as fraction of total cropland

# import data
raa <- "InData/CropAbundance/*.tif" %>% Sys.glob() %>% rast()
raa <- classify(raa, cbind(0, NA))
crops <- names(raa)

tcl <- rast("InData/TotalCropland.tif")

dir <- "Maps/ActualCropArea"
dir.create(dir, F, T)
ar = ncol(raa[[1]])/nrow(raa[[1]])

for(i in 1:length(crops)){
  r <- lapp(c(raa[[i]], tcl), `/`)
  brks <- c(0.0001, 0.001, 0.01, 0.1, 0.2, 0.4, 0.8, 1)
  fgfn = file.path(dir, paste0(crops[i], '.tif'))
  tiff(filename = fgfn, units = "in",
       width = ncol(abund[[i]])/300, 
       height = (ncol(abund[[i]])/300)/ar, 
       type = "cairo", res = 300, 
       compression = "zip")
  plot(r, breaks = brks, col = RColorBrewer::brewer.pal(length(brks), "Greens"),
       axes = FALSE, maxcell = ncell(r), 
       main = paste(crops[i], "actural area (fraction"))
  plot(countries, lwd = 1, add = T, border = 'red', lwd = 0.2)
  dev.off()
}

# Suitabilities ################################################################

## Maxent ----
fn <- Sys.glob("OutData/SDM/Suitability/*.tif")
me_suit <- rast(fn)
me_dir <- "Maps/Suitability/Maxent"
dir.create(me_dir, F, T)
crops <- gsub("OutData/SDM/Suitability/", "",
              gsub(".suit.tif", "", fn))

ar = ncol(me_suit[[1]])/nrow(me_suit[[1]])

for(i in 1:length(crops)){
  { 
    fgfn = file.path(me_dir, paste0(crops[i], 'ME.tif'))
    tiff(filename = fgfn, units = "in",
         width = ncol(me_suit[[i]])/300, 
         height = (ncol(me_suit[[i]])/300)/ar, 
         type = "cairo", res = 300, 
         compression = "zip")
    {
      plot(me_suit[[i]], axes = FALSE, maxcell = ncell(me_suit[[i]]), 
           main = crops[i])
      plot(ne_countries(), lwd = 1, add = T)
    }
    dev.off()
  }
}


## Ecocrop ----

crops <- c('maize', 'wheat', 'rice', 'soybean', 'oilpalm', 
           'avocado', 'grape', 'apple', 'almond')

sfn <- Sys.glob("D:/globcropdiv/EcocropSuit/byMonfCat/*.tif")
suit <- rast(sfn)
smdir <- "Maps/Suitability"
dir.create(smdir)
smdir <- "Maps/Suitability/EcoCrop"
dir.create(smdir)


ar = ncol(suit[[1]])/nrow(suit[[1]])

for(i in 1:length(crops)){
  { # run this line to save plot 
    ci <- which(names(suit) == crops[i])
    fgfn = file.path(smdir, paste0(crops[i], '.tif'))
    tiff(filename = fgfn, units = "in",
         width = ncol(suit[[ci]])/300, 
         height = (ncol(suit[[ci]])/300)/ar, 
         type = "cairo", res = 300, 
         compression = "zip")
    {
      plot(suit[[ci]], axes = FALSE, 
           maxcell = ncell(suit[[ci]]), main = crops[i])
      plot(ne_countries(), lwd = 1, add = T)
    }
    dev.off()
  }
}


## Integrated ----
fn <- Sys.glob("OutData/Int_Suit/*.tif") 
suit <- rast(fn)
dir <- "Maps/Suitability/Integrated"
dir.create(dir, F, T)
crops <- gsub("OutData/Int_Suit/", "",
              gsub(".suit.tif", "", fn))

ar = ncol(suit[[1]])/nrow(suit[[1]])

for(i in 1:length(crops)){
  { 
    fgfn = file.path(dir, paste0(crops[i], '_suit.tif'))
    tiff(filename = fgfn, units = "in",
         width = ncol(suit[[i]])/300, 
         height = (ncol(suit[[i]])/300)/ar, 
         type = "cairo", res = 300, 
         compression = "zip")
    {
      plot(suit[[i]], axes = FALSE, maxcell = ncell(suit[[i]]), 
           main = crops[i])
      plot(ne_countries(), lwd = 1, add = T)
    }
    dev.off()
  }
}
