library(terra)

fn <- Sys.glob("D:/Monfreda/GeoTiff/*/*HarvestedAreaHectares.tif")
monf_area <- rast(fn)

f <- function(x){
  if(all(is.na(x) | x == 0)) return(NA)
  pi <- x/(sum(x, na.rm = T))
  D <- exp(-sum(log(pi) * pi, na.rm = T))
  return(D)
}

app(monf_area, fun = f, filename = "D:/Dp_World/LowRes_10min/Da_Monf_5min.tif",
    nodes = 7, overwrite = T, 
    wopt = list(names = "Da_Monf.tif", filetype = "GTiff", progress = 1))

########### 
Da <- rast("D:/Dp_World/LowRes_10min/Da_Monf_5min.tif")

plot(Da, col = colorRamps::matlab.like2(400))
plot(rnaturalearth::ne_countries(), add  = T)
