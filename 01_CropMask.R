library(terra)

nwd <- nchar(getwd())
if(substr(getwd(),  nwd - 10, nwd) != "globcropdiv") stop("See 0000_wd.R")

# upload all Monfreda Crop Proportion Data
Ddrive <- dir.exists("D:/Monfreda") 
if(Ddrive) {
  fn <- Sys.glob("D:/Monfreda/GeoTiff/*/*HarvestedAreaFraction.tif")
} else {
  fn <- Sys.glob("InData/Monfreda/GeoTiff/*/*HarvestedAreaFraction.tif")
}

# exclude Crops that aren't in Ecocrops (see 00_CropSpeciesSelection.R)
torm <- read.csv("AuxData/ExcMonfCrops.csv")
for(i in 1:nrow(torm)){
  fn <- fn[!grepl(torm$Monf_Name[i], fn)]  
}

monf_area <- rast(fn)

# create cropland mask
f <- function(x) ifelse(all(x == 0), NA, 1)

if(Ddrive) dir.create("D:/globcropdiv") else dir.create("OutData")

outfn <- file.path(ifelse(Ddrive, "D:/globcropdiv", "OutData"),
                   "Monf_Cropland_Mask.tif")

app(monf_area, fun = f, 
    filename = outfn,
    nodes = 7, overwrite = T, 
    wopt = list(names = "CropMask", filetype = "GTiff", progress = 1))

cm <- rast(outfn)
plot(cm)
