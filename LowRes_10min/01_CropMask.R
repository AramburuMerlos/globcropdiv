library(terra)
library(Recocrop)

# weather data
tavg <- rast(Sys.glob("D:/WorldClim/2.1/wc10min/tavg/*.tif")) 

# function to aggregate cropland mask to same resolution as weather data
# non-cropland (including water) only when there is no cropland (value 0)
# otherwise, majority count

# cropland mask
clm <- rast("D:/GlobalCropMask_NASA/GFSAD1KCM.2010.001.2016348142550.tif")
# 9 is non-cropland areas
# 1 is (major) irrigated cropland areas. 
# 2 is (minor) irrigated cropland areas.
# Maybe for 1 and 2 I can simulate crop suitabilities without accounting for rainfall 

# In the meantime I'll consider all cropland areas

f1 <- function(x) ifelse(all(x == 0 | x == 9), 0, which.max(tabulate(x)[1:5]))
agfact <- floor(res(tavg) / res(clm))
clm <- terra::aggregate(clm, fact = agfact, fun = f1, nodes = 8)
resample(clm, tavg[[1]], method = "near", overwrite = T,
         filename = "D:/GlobalCropMask_NASA/CropMask10min.tif")

