library(raster)
library(terra)
library(data.table)

# Dp ############
fn <- Sys.glob("D:/Dp_world/CropSuit/byMonfCat/*.tif")
suit <- rast(fn)

f <- function(x){
  if(all(is.na(x))) return(NA)
  pi <- x/(sum(x, na.rm = T))
  D <- exp(-sum(log(pi) * pi, na.rm = T))
  return(D)
}

app(suit, fun = f, filename = "D:/Dp_World/Dp_MonfCat.tif",
    nodes = 7, overwrite = T, 
    wopt = list(names = "Dp_Monf", filetype = "GTiff", progress = 1))

# Da ############
fn <- Sys.glob("D:/Monfreda/GeoTiff/*/*HarvestedAreaFraction.tif")
monf_area <- rast(fn)

f <- function(x){
  if(all(is.na(x) | x == 0)) return(NA)
  pi <- x/(sum(x, na.rm = T))
  D <- exp(-sum(log(pi) * pi, na.rm = T))
  return(D)
}

app(monf_area, fun = f, filename = "D:/Dp_World/Da_MonfCat.tif",
    nodes = 7, overwrite = T, 
    wopt = list(names = "Da_Monf.tif", filetype = "GTiff", progress = 1))

# Maps #########
Dp <- rast("D:/Dp_World/Dp_MonfCat.tif")
Da <- rast("D:/Dp_World/Da_MonfCat.tif")

plot(Dp, col = colorRamps::matlab.like2(400))
plot(rnaturalearth::ne_countries(), add  = T)


plot(Da, col = colorRamps::matlab.like2(400))
plot(rnaturalearth::ne_countries(), add  = T)

Dpf <- app(c(Da,Dp), max, na.rm = T)

plot(Da/Dpf, col = colorRamps::matlab.like2(400))
plot(rnaturalearth::ne_countries(), add  = T)



Dpv <- getValues(Dp)
Dav <- getValues(Da)

DT <- data.table(Dp = Dpv, Da = Dav)

DT <- DT[!is.na(Dp) & !is.na(Da),]
DT[, Ha:= log(Da)]
DT[, Hp:= log(Dp)]

plot(DT$Ha, DT$Hp, xlim = c(0,6), ylim = c(0,6))
DT[, cor(Ha, Hp)]
