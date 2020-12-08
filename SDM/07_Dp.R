library(raster)
library(terra)
library(data.table)
library(rnaturalearth)

nwd <- nchar(getwd())
if(substr(getwd(),  nwd - 10, nwd) != "globcropdiv") warning("See 0000_wd.R")

path <- "OutData/SDM"


###################        Maxent       #######################################
# Dp ----------
fn <- Sys.glob(file.path(path, "Maxent/*.tif"))
suit <- rast(fn)
crops <- sub(paste0(path, "/Maxent/"), "", sub("_Maxent_pred.tif", "", fn))

f <- function(x){
  x[x==0] <- NA
  if(all(is.na(x))) return(NA)
  pi <- x/(sum(x, na.rm = T))
  D <- exp(-sum(log(pi) * pi, na.rm = T))
  return(D)
}

app(suit, fun = f, filename = "OutData/Dp_Maxent.tif",
    overwrite = T, wopt = list(names = "Dp_Maxent", filetype = "GTiff", progress = 1))


# Map plot --------
Dp <- rast("OutData/Dp_Maxent.tif")

ar = ncol(Dp)/nrow(Dp)

{ # run this line to save plot 
  fgfn = ('Maps/Dp_Maxent.tif')
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
         mar = c(2,0,2,4), main = "'Potential' Diversity (Maxent)")
    plot(ne_countries(), lwd = 1, add = T)
  }
  dev.off()
}

####### Suitabilities ################
sc <- fread("AuxData/CropAbundanceSource_Filtered.csv")

smdir <- "Maps/Suitability/Maxent"
dir.create(smdir)

ar = ncol(suit[[1]])/nrow(suit[[1]])

for(i in 1:length(crops)){
  { # run this line to save plot 
    fgfn = file.path(smdir, paste0(crops[i], '.tif'))
    tiff(filename = fgfn, units = "in",
         width = ncol(suit[[i]])/300, 
         height = (ncol(suit[[i]])/300)/ar, 
         type = "cairo", res = 300, 
         compression = "zip")
    {
      plot(suit[[i]], axes = FALSE, maxcell = ncell(suit[[i]]), 
           main = sc$FAO_Name[sc$filename == crops[i]])
      plot(ne_countries(), lwd = 1, add = T)
    }
    dev.off()
  }
}

