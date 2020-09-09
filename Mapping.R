library(terra)
library(rnaturalearth)

Dp <- rast("D:/Dp_World/Dp_MonfCat.tif")
Da <- rast("D:/Dp_World/Da_MonfCat.tif")

dir.create("Maps")

ar = ncol(Dp)/nrow(Dp)

{ # run this line to save plot 
  fgfn = ('Maps/Dp_MonfCat.tif')
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
         mar = c(2,0,2,4))
    plot(ne_countries(), lwd = 1, add = T)
  }
  dev.off()
}

{ # run this line to save plot 
  fgfn = ('Maps/Da_MonfCat.tif')
  tiff(filename = fgfn, units = "in",
       width = ncol(Da)/300, 
       height = (ncol(Da)/300)/ar, 
       type = "cairo", res = 300, 
       compression = "zip")
  {
    par(mai = c(0,0,0,0), omi = c(0,0,0,0))
    cols = colorRamps::matlab.like2(400)
    plot(Da, axes = FALSE, col = cols,
         maxcell = ncell(Da), 
         mar = c(2,0,2,4))
    plot(ne_countries(), lwd = 1, add = T)
  }
  dev.off()
}

