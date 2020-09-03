library(terra)

fn <- Sys.glob("D:/Dp_world/LowRes_10min/CropSuit/*.tif")
suit <- rast(fn)

f <- function(x){
  if(all(is.na(x))) return(NA)
  pi <- x/(sum(x, na.rm = T))
  D <- exp(-sum(log(pi) * pi, na.rm = T))
  return(D)
}

app(suit, fun = f, filename = "D:/Dp_World/LowRes_10min/Dp_10min.tif",
    nodes = 7, overwrite = T, 
    wopt = list(names = "Dp_10min", filetype = "GTiff", progress = 1))

########### 
Dp <- rast("D:/Dp_World/LowRes_10min/Dp_10min.tif")

plot(Dp, col = colorRamps::matlab.like2(400))
plot(rnaturalearth::ne_countries(), add  = T)
summary(Dp) 
