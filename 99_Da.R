library(raster)
library(terra)
library(data.table)
library(rnaturalearth)


nwd <- nchar(getwd())
if(substr(getwd(),  nwd - 10, nwd) != "globcropdiv") warning("See 0000_wd.R")

# directory of input data
indir <- "D:"

countries <- ne_download(scale = 10, type = "countries")

# Crop Abundance Data ----------
sc <- fread("AuxData/CropAbundanceSource.csv")
crops <- sc[,ifelse(source == "SPAM", SPAM_Name, Monf_Name)]
abund <- rast(file.path(indir, sc$file))
names(abund) <- crops

####### Abundance Maps ################
dir <- "Maps/Abundance"
dir.create(dir)

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


######## Actual Diversity ##########
f <- function(x){
  if(all(is.na(x) | x == 0)) return(NA)
  pi <- x/(sum(x, na.rm = T))
  D <- exp(-sum(log(pi) * pi, na.rm = T))
  return(D)
}

app(abund, fun = f, filename = "OutData/Da_all.tif",
    overwrite = T, 
    wopt = list(names = "Da", filetype = "GTiff", progress = 1))

# Map plot --------
Da <- rast("OutData/Da_all.tif")
ar = ncol(Da)/nrow(Da)

{ # run this line to save plot 
  fgfn = ('Maps/Da_all.tif')
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
         mar = c(2,0,2,4), main = "'Actual Diversity (all crops)")
    plot(countries, lwd = 0.4, add = T)
  }
  dev.off()
}
