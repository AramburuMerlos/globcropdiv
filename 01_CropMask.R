library(terra)
# upload all Monfreda Crop Proportion Data
fn <- Sys.glob("D:/Monfreda/GeoTiff/*/*HarvestedAreaFraction.tif")

# exclude Crops that aren't in Ecocrops (see 00_CropSpeciesSelection.R)
torm <- read.csv("AuxData/ExcMonfCrops.csv")
for(i in 1:nrow(torm)){
  fn <- fn[!grepl(torm$Monf_Name[i], fn)]  
}

monf_area <- rast(fn)

# create cropland mask
f <- function(x) ifelse(all(x == 0), NA, 1)

app(monf_area, fun = f, filename = "D:/Dp_World/Monf_Cropland_Mask.tif",
    nodes = 7, overwrite = T, 
    wopt = list(names = "CropMask", filetype = "GTiff", progress = 1))

cm <- rast("D:/Dp_World/Monf_Cropland_Mask.tif")
plot(cm)
