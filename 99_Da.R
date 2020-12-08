library(raster)
library(terra)
library(data.table)
library(rnaturalearth)


nwd <- nchar(getwd())
if(substr(getwd(),  nwd - 10, nwd) != "globcropdiv") warning("See 0000_wd.R")

# directory of input data
indir <- "D:"
cpath <- file.path(indir, "WorldCropAbundance")
path <- "OutData/SDM"

# Crop Abundance Data ----------
sc <- fread("AuxData/CropAbundanceSource_Filtered.csv")
# KEEP ONLY CROPS WITH >75% Of THEIR DATA NON-FILTERED 
sel_sc <- sc[nf_prop_ha > 0.75,]
crops <- sel_sc[, filename]
crops.fn <- paste0(crops, "_", sel_sc[, source])
fn <- file.path(cpath, paste0("fraction/", crops.fn, "_fr.tif"))
abund <- rast(fn)

countries <- ne_download(scale = 10, type = "countries")
####### Abundance Maps ################
dir <- "Maps/Abundance/Filtered"
dir.create(dir)

ar = ncol(abund[[1]])/nrow(abund[[1]])

for(i in 1:length(crops)){
  { # run this line to save plot 
    fgfn = file.path(dir, paste0(crops[i], '.tif'))
    tiff(filename = fgfn, units = "in",
         width = ncol(abund[[i]])/300, 
         height = (ncol(abund[[i]])/300)/ar, 
         type = "cairo", res = 300, 
         compression = "zip")
    {
      plot(raster(abund[[i]]), 
           axes = FALSE, maxpixels = ncell(abund[[i]]), main = sel_sc$FAO_Name[i], colNA = 'black')
      plot(countries, lwd = 1, add = T, border = 'red', lwd = 0.2)
    }
    dev.off()
  }
}

######## Actual Diversity ##########
#with filtered Data

f <- function(x){
  if(all(is.na(x) | x == 0)) return(NA)
  pi <- x/(sum(x, na.rm = T))
  D <- exp(-sum(log(pi) * pi, na.rm = T))
  return(D)
}

app(abund, fun = f, filename = "OutData/Da_88crops_filtered.tif",
    overwrite = T, 
    wopt = list(names = "Da", filetype = "GTiff", progress = 1))

# Map plot --------
Da <- rast("OutData/Da_88crops_filtered.tif")
ar = ncol(Da)/nrow(Da)

{ # run this line to save plot 
  fgfn = ('Maps/Da_88crops_filtered.tif')
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
         mar = c(2,0,2,4), main = "'Actual Diversity (88 crops)")
    plot(countries, lwd = 0.4, add = T)
  }
  dev.off()
}
