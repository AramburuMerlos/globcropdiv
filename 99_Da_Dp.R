library(magrittr)
library(raster)
library(terra)
library(data.table)
library(rnaturalearth)

if(!getwd() %like%  "globcropdiv$") warning("See 0000_wd.R")

countries <- ne_download(scale = 10, type = "countries")

f <- function(x){
  if(all(is.na(x) | x == 0)) return(NA)
  pi <- x/(sum(x, na.rm = T))
  D <- exp(-sum(log(pi) * pi, na.rm = T))
  return(D)
}

# cells with cropland
totcl <- rast("InData/TotalCropland.tiff")
cells <- values(totcl) %>% is.na %>% `!` %>% which()


# Actual Diversity #############################################################

# Crop Abundance Data
Da <- "AuxData/CropAbundanceSource.csv" %>% 
  fread %$% file %>% file.path("D:", .) %>% 
  rast %>% extract(cells) %>% apply(1,f)

# assign values
r <- rast(totcl)
v <- rep(NA, ncell(r))
v[cells] <- Da
values(r) <- v

writeRaster(r, filename = "OutData/Da.tif",
            overwrite = T, 
            wopt = list(names = "Da", filetype = "GTiff", progress = 1))


# Potential Diversity ##########################################################

# Maxent suitability -------------
Dp <- "OutData/SDM/Suitability/*.tif" %>%
  Sys.glob() %>% rast %>% extract(cells) %>% apply(1, f)

# assign values
r <- rast(totcl)
v <- rep(NA, ncell(r))
v[cells] <- Dp
values(r) <- v

writeRaster(r, filename = "OutData/Dp_MEsuit.tif",
            overwrite = T, 
            wopt = list(names = "Dp_MEsuit", filetype = "GTiff", progress = 1))


# Ecocrop Suitability ----------------------
Dp <- "OutData/EcocropSuit/byMonfCat/*.tif" %>%
  Sys.glob %>% rast %>% extract(cells) %>% apply(1,f)


# assign values
r <- rast(totcl)
v <- rep(NA, ncell(r))
v[cells] <- Dp
values(r) <- v

writeRaster(r, filename = "OutData/Dp_Ecocrop.tif",
            overwrite = T, 
            wopt = list(names = "Dp_MEsuit", filetype = "GTiff", progress = 1))

