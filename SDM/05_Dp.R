library(raster)
library(terra)
library(data.table)
library(rnaturalearth)

nwd <- nchar(getwd())
if(substr(getwd(),  nwd - 10, nwd) != "globcropdiv") stop("See 0000_wd.R")

path <- "D:/globcropdiv"
path <- ifelse(dir.exists(path), path, "OutData") 


# Dp ############
fn <- Sys.glob(file.path(path, "SDM/*.tif"))
suit <- rast(fn)
crops <- sub(paste0(path, "/SDM/"), "", sub("_RFpred.tif", "", fn))

f <- function(x){
  x[x==0] <- NA
  if(all(is.na(x))) return(NA)
  pi <- x/(sum(x, na.rm = T))
  D <- exp(-sum(log(pi) * pi, na.rm = T))
  return(D)
}

app(suit, fun = f, filename = file.path(path, "Dp_RF_MonfCat.tif"),
    nodes = 7, overwrite = T, 
    wopt = list(names = "Dp_Monf", filetype = "GTiff", progress = 1))


# Map plot --------
Dp <- rast(file.path(path, "Dp_RF_MonfCat.tif"))

ar = ncol(Dp)/nrow(Dp)

{ # run this line to save plot 
  fgfn = ('Maps/Dp_RF_MonfCat.tif')
  tiff(filename = fgfn, units = "in",
       width = ncol(Dp)/300, 
       height = (ncol(Dp)/300)/ar, 
       type = "cairo", res = 300, 
       compression = "zip")
  {
    par(mai = c(0,0,0,0), omi = c(0,0,0,0))
    cols = colorRamps::matlab.like2(400)
    plot(Dp, axes = FALSE, col = cols,
         maxcell = ncell(Dp), 
         mar = c(2,0,2,4), main = "'Potential' Diversity (RF - Monfreda)")
    plot(ne_countries(), lwd = 1, add = T)
  }
  dev.off()
}

####### Suitabilities ################
crops <- c('maize', 'wheat', 'rice', 'soybean', 'oilpalm', 
           'avocado', 'grape', 'apple', 'almond', 'clover')

smdir <- "Maps/Suitability/RF"
dir.create(smdir)

ar = ncol(suit[[1]])/nrow(suit[[1]])

for(i in 1:length(crops)){
  { # run this line to save plot 
    ci <- which(names(suit) == crops[i])
    fgfn = file.path(smdir, paste0(crops[i], '.tif'))
    tiff(filename = fgfn, units = "in",
         width = ncol(suit[[ci]])/300, 
         height = (ncol(suit[[ci]])/300)/ar, 
         type = "cairo", res = 300, 
         compression = "zip")
    {
      plot(suit[[ci]], axes = FALSE, maxcell = ncell(suit[[ci]]), main = crops[i])
      plot(ne_countries(), lwd = 1, add = T)
    }
    dev.off()
  }
}
