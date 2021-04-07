library(data.table)
library(raster)
library(terra)
library(rnaturalearth)
library(colorRamps)

if(!getwd() %like%  "globcropdiv$") warning("See 0000_wd.R")

countries <- ne_download(scale = 10, type = "countries")

# Irrigation ###################################################################
# check FAO Aquastat data
iperc <- rast("D:/AQUASTAT/gmia_v5_aei_pct.asc")

iperc <- classify(iperc, rcl = cbind(0,NA))

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
             'cadetblue2', 'dodgerblue2', 'mediumpurple3', 
             'pink2', 'firebrick3')
    plot(iperc, axes = FALSE, type = "interval", maxcell = ncell(iperc),
         col = cols, levels = c(0,0.1,1,5,10,20,35,50,75,100),
         mar = c(2,0,3,5), main = "Area Equiped for Irrigation")
    plot(ne_countries(), lwd = 1, add = T)
  }
  dev.off()
}

# Abundance Maps ###############################################################

# import data
sc <- fread("AuxData/CropAbundanceSource.csv")
crops <- sc[,ifelse(source == "SPAM", SPAM_Name, Monf_Name)]
abund <- rast(file.path("D:", sc$file))
names(abund) <- crops

dir <- "Maps/Abundance"
dir.create(dir, F, T)
ar = ncol(abund[[1]])/nrow(abund[[1]])

for(i in 1:length(crops)){
  r <- raster(abund[[i]])
  r <- reclassify(r, rcl = cbind(0,NA))
  rmax <- cellStats(r, 'max')
  brks <- c(0.01,1,if(rmax > 1000) c(100, 1000, rmax) 
            else if (rmax>100) c(100, rmax) else rmax)
  fgfn = file.path(dir, paste0(crops[i], '.tif'))
  tiff(filename = fgfn, units = "in",
       width = ncol(abund[[i]])/300, 
       height = (ncol(abund[[i]])/300)/ar, 
       type = "cairo", res = 300, 
       compression = "zip")
  plot(r, breaks = brks, col = RColorBrewer::brewer.pal(length(brks), "Greens"),
       axes = FALSE, maxpixels = ncell(r), 
       main = sc$FAO_Name[i], colNA = 'black')
  plot(countries, lwd = 1, add = T, border = 'red', lwd = 0.2)
  dev.off()
}


# data quality index for Monfreda crops
dqi.fn <- Sys.glob("D:/Monfreda/DataQualityIndex/*.tif")
crops_dqi <- gsub("D:/Monfreda/DataQualityIndex/", "", 
                  gsub(".tif", "", dqi.fn))
dqi <- rast(dqi.fn)

ar = ncol(dqi[[1]])/nrow(dqi[[1]])

for(i in 1:length(crops_dqi)){
  fgfn = file.path(dir, paste0(crops_dqi[i], '.tif'))
  tiff(filename = fgfn, units = "in",
       width = ncol(dqi[[i]])/300, 
       height = (ncol(dqi[[i]])/300)/ar, 
       type = "cairo", res = 300, 
       compression = "zip")
  plot(dqi[[i]], col = colorRamps::green2red(100), maxcell = ncell(dqi[[i]]), 
       main = sc[source=="Monf", FAO_Name][i])
  plot(countries, lwd = 1, add = T, lwd = .5)
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


# Actual and Potential Diversity ##############################################
Da <- rast("OutData/Da.tif")
Dp_me <- rast("OutData/Dp_MEsuit.tif")
Dp_ec <- rast("OutData/Dp_Ecocrop.tif")

ar = ncol(Da)/nrow(Da)

{ # run this line to save plot 
  fgfn = ('Maps/Da.tif')
  tiff(filename = fgfn, units = "in",
       width = ncol(Da)/300, 
       height = (ncol(Da)/300)/ar, 
       type = "cairo", res = 300, 
       compression = "zip")
  {
    par(mai = c(0,0,0,0), omi = c(0,0,0,0))
    brks = c(seq(1,10,1), seq(15,50,5), seq(60,170,10))
    cols = rev(matlab.like2(length(brks)-1))
    plot(Da, type = "interval", axes = FALSE, col = cols, breaks = brks,
         maxcell = ncell(Da), 
         mar = c(2,0,2,4), main = "Actual Diversity")
    plot(countries, lwd = 0.4, add = T)
  }
  dev.off()
}

{ # run this line to save plot 
  fgfn = ('Maps/Dp_ecocrop.tif')
  tiff(filename = fgfn, units = "in",
       width = ncol(Dp_ec)/300, 
       height = (ncol(Dp_ec)/300)/ar, 
       type = "cairo", res = 300, 
       compression = "zip")
  {
    par(mai = c(0,0,0,0), omi = c(0,0,0,0))
    brks = c(seq(1,10,1), seq(15,50,5), seq(60,170,10))
    cols = rev(matlab.like2(length(brks)-1))
    plot(Dp_ec, type = "interval", axes = FALSE, col = cols, breaks = brks,
         maxcell = ncell(Dp_ec), 
         mar = c(2,0,2,4), main = "Potential Diversity (Ecocrop)")
    plot(countries, lwd = 0.4, add = T)
  }
  dev.off()
}

{ # run this line to save plot 
  fgfn = ('Maps/Dp_MEsuit.tif')
  tiff(filename = fgfn, units = "in",
       width = ncol(Dp_me)/300, 
       height = (ncol(Dp_me)/300)/ar, 
       type = "cairo", res = 300, 
       compression = "zip")
  {
    par(mai = c(0,0,0,0), omi = c(0,0,0,0))
    brks = c(seq(1,10,1), seq(15,50,5), seq(60,170,10))
    cols = rev(matlab.like2(length(brks)-1))
    plot(Dp_me, type = "interval", axes = FALSE, col = cols, breaks = brks,
         maxcell = ncell(Dp_me), 
         mar = c(2,0,2,4), main = "Potential Diversity (ME suit)")
    plot(countries, lwd = 0.4, add = T)
  }
  dev.off()
}


