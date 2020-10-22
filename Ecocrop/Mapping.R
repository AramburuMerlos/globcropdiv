library(terra)
library(rnaturalearth)

setwd(here::here())
nwd <- nchar(getwd())
if(substr(getwd(),  nwd - 10, nwd) != "globcropdiv") stop("See 0000_wd.R")


#### Diversity (Monfreda Categ) ###############
Dp <- rast("D:/globcropdiv/Dp_MonfCat.tif")
Da <- rast("D:/globcropdiv/Da_MonfCat.tif")

dir.create("Maps")

ar = ncol(Dp)/nrow(Dp)

{ # run this line to save plot 
  fgfn = ('Maps/Dp_EcoCrop_MonfCat.tif')
  tiff(filename = fgfn, units = "in",
       width = ncol(Dp)/300, 
       height = (ncol(Dp)/300)/ar, 
       type = "cairo", res = 300, 
       compression = "zip")
  {
    par(mai = c(0,0,0,0), omi = c(0,0,0,0))
    cols = colorRamps::matlab.like2(400)
    plot(Dp, axes = FALSE, col = cols,
         maxcell = ncell(Dp), 
         mar = c(2,0,2,4), main = "'Potential' Diversity (Monfreda)")
    plot(ne_countries(), lwd = 1, add = T)
  }
  dev.off()
}

{ # run this line to save plot 
  fgfn = ('Maps/Da_MonfCat.tif')
  tiff(filename = fgfn, units = "in",
       width = ncol(Da)/300, 
       height = (ncol(Da)/300)/ar, 
       type = "cairo", res = 300, 
       compression = "zip")
  {
    par(mai = c(0,0,0,0), omi = c(0,0,0,0))
    cols = colorRamps::matlab.like2(400)
    plot(Da, axes = FALSE, col = cols,
         maxcell = ncell(Da), 
         mar = c(2,0,2,4), main = "Actual Diversity (Monfreda)")
    plot(ne_countries(), lwd = 1, add = T)
  }
  dev.off()
}

####### Irrigation ##############
# check FAO Aquastat data
iperc <- rast("D:/AQUASTAT/gmia_v5_aei_pct.asc")
iperc <- app(iperc, fun = function(x)ifelse(x == 0, NA, x), nodes = 5)
ar = ncol(iperc)/nrow(iperc)

{ # run this line to save plot 
  fgfn = ('Maps/FAOP_AQUASTAT.tif')
  tiff(filename = fgfn, units = "in",
       width = ncol(iperc)/300, 
       height = (ncol(iperc)/300)/ar, 
       type = "cairo", res = 300, 
       compression = "zip")
  {
    par(mai = c(0,0,0,0), omi = c(0,0,0,0))
    cols = c('khaki1', 'goldenrod2', 'palegreen', 'darkgreen',
             'cadetblue2', 'dodgerblue2', 'mediumpurple3', 'pink2', 'firebrick3')
    plot(iperc, axes = FALSE, type = "interval", maxcell = ncell(iperc),
         col = cols, levels = c(0,0.1,1,5,10,20,35,50,75,100),
         mar = c(2,0,3,5), main = "Area Equiped for Irrigation")
    plot(ne_countries(), lwd = 1, add = T)
  }
  dev.off()
}

####### Suitabilities ################
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
      plot(suit[[ci]], axes = FALSE, maxcell = ncell(suit[[ci]]), main = crops[i])
      plot(ne_countries(), lwd = 1, add = T)
    }
    dev.off()
  }
}
  
####### Areas ################
crops <- c('maize', 'wheat', 'rice', 'soybean', 'oilpalm', 
           'avocado', 'grape', 'apple', 'almond')

cmask <- rast("OutData/Monf_Cropland_Mask.tif")
afn <- Sys.glob("D:/Monfreda/GeoTiff/*/*HarvestedAreaFraction.tif")
monf_area <- rast(afn)
monf_names <- sub("_HarvestedAreaFraction", "", names(monf_area))

adir <- "Maps/Area"
dir.create(adir)

ar = ncol(monf_area[[1]])/nrow(monf_area[[1]])

for(i in 1:length(crops)){
  { # run this line to save plot 
    ci <- which(crops[i] == monf_names)
    fgfn = file.path(adir, paste0(crops[i], '.tif'))
    tiff(filename = fgfn, units = "in",
         width = ncol(monf_area[[ci]])/300, 
         height = (ncol(monf_area[[ci]])/300)/ar, 
         type = "cairo", res = 300, 
         compression = "zip")
    {
      r <- mask(monf_area[[ci]], cmask)
      v <- getValues(raster(r))
      v <- v[!is.na(v) & v>0]
      brks <- c(0, quantile(v,c(0.5,0.75,1)))
      plot(r, axes = FALSE, maxcell = ncell(monf_area[[ci]]), 
           main = crops[i], type = "interval", levels = brks, 
           col = rev(terrain.colors(length(brks)-1)))
      plot(ne_countries(), lwd = 1, add = T)
    }
    dev.off()
  }
}



