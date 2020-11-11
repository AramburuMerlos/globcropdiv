library(data.table)
library(raster)
library(terra)
library(rnaturalearth)


setwd(here::here())
nwd <- nchar(getwd())
if(substr(getwd(),  nwd - 10, nwd) != "globcropdiv") stop("See 0000_wd.R")

# directory or drive where the raw data is storage
rd <- "D:/"

####### Crop Abundance ################
# for selected crops
sc <- fread("AuxData/CropAbundanceSource.csv")
# SPAM ----
# upload SPAM data
scode <- sc[source == "SPAM", SPAM_Code]
sfn <- vector(length = length(scode))
for(i in 1:length(scode)){
  sfn[i] <- paste0(rd, "SPAM/Physical_Area/spam2010V2r0_global_A_", 
                   toupper(scode[i]), "_A.tif")
}
scrops <- rast(sfn)
names(scrops) <- paste0("SPAM_", sc[source == "SPAM", SPAM_Name], "_(ha)")

# change 0 to NA
scrops <- classify(scrops, cbind(0, NA), othersNA = FALSE)


cadir <- "Maps/CropAbundance"
dir.create(cadir)

ar = ncol(scrops[[1]])/nrow(scrops[[1]])

for(i in 1:nlyr(scrops)){
  { # run this line to save plot
    fgfn = file.path(cadir, paste0(names(scrops[[i]]), '.tif'))
    tiff(filename = fgfn, units = "in",
         width = ncol(scrops[[i]])/300, 
         height = (ncol(scrops[[i]])/300)/ar, 
         type = "cairo", res = 300, 
         compression = "zip")
    {
      v <- getValues(raster(scrops[[i]]))
      v <- v[!is.na(v)]
      brks <- c(0, quantile(v,c(0.5,0.75,1)))
      plot(scrops[[i]], axes = FALSE, maxcell = ncell(scrops[[i]]), 
           main = names(scrops[[i]]), type = "interval", levels = brks, 
           col = rev(terrain.colors(length(brks)))[2:length(brks)])
      plot(ne_countries(), lwd = 1, add = T)
    }
    dev.off()
  }
}







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
  
