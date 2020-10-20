library(raster)
library(terra)
library(data.table)

nwd <- nchar(getwd())
if(substr(getwd(),  nwd - 10, nwd) != "globcropdiv") stop("See 0000_wd.R")

path <- "D:/globcropdiv"
path <- ifelse(dir.exists(path), path, "OutData") 


# Dp ############
fn <- Sys.glob(file.path(path, "EcocropSuit/byMonfCat/*.tif"))
suit <- rast(fn)

f <- function(x){
  if(all(is.na(x))) return(NA)
  pi <- x/(sum(x, na.rm = T))
  D <- exp(-sum(log(pi) * pi, na.rm = T))
  return(D)
}

app(suit, fun = f, filename = file.path(path, "Dp_EcoCrop_MonfCat.tif"),
    nodes = 7, overwrite = T, 
    wopt = list(names = "Dp_Monf", filetype = "GTiff", progress = 1))

# Da ############
monpath <- if(dir.exists("D:/Monfreda")) "D:/Monfreda" else ("InData/Monfreda")
fn <- Sys.glob(file.path(monpath, "/GeoTiff/*/*HarvestedAreaFraction.tif"))
monf_area <- rast(fn)

f <- function(x){
  if(all(is.na(x) | x == 0)) return(NA)
  pi <- x/(sum(x, na.rm = T))
  D <- exp(-sum(log(pi) * pi, na.rm = T))
  return(D)
}

app(monf_area, fun = f, filename = file.path(path, "Da_MonfCat.tif"),
    nodes = 7, overwrite = T, 
    wopt = list(names = "Da_Monf.tif", filetype = "GTiff", progress = 1))

