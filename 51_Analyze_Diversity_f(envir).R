library(data.table)
library(terra)


if(system('hostname', TRUE) == "ESP-RH-9891"){
  setwd("D:/globcropdiv/")
} else if(system('hostname', TRUE) == "LAPTOP-ST129J47"){ 
  setwd("G:/My Drive/globcropdiv/")
} # else if { ... 

# upload data --------------

totcl <- rast("InData/TotalCropland.tif")
area <- cellSize(totcl, unit = "ha")
pcl <- totcl/area
cl_mask <- pcl < 0.005
totcl <- mask(totcl, cl_mask, maskvalue = 1)


rGDD <- rast("InData/WorldClim/2.1/wc5min/extra/GDD.tif")
rtm <- rast("InData/WorldClim/2.1/wc5min/bioc/wc2.1_5m_bio_1.tif")
rpp <- rast("InData/WorldClim/2.1/wc5min/bioc/wc2.1_5m_bio_12.tif")
rai <- rast("InData/WorldClim/2.1/wc5min/extra/AI.tif")


# extract values
d <- data.table(cells = which(!is.na(values(totcl))))
d[, tcl:= extract(totcl, cells)[,1]]
d[, GDD:= extract(rGDD, cells)[,1]/365]
d[, tm:= extract(rtm, cells)[,1]]
d[, pp:= extract(rpp, cells)[,1]]
d[, ai:= extract(rai, cells)[,1]]
d[ai > 3, ai:= 3]

rdivs <- rast(Sys.glob("OutData/*_D*tif"))
divs <- names(rdivs)
d[, (divs):= extract(rdivs, cells)]


# create breaks ---------------------
nbreaks <- 51
fbrks <- function(x) seq(min(x, na.rm = T), max(x, na.rm = T), l = nbreaks)
dbrks <- d[, lapply(.SD, fbrks) , .SDcols = divs]
gdd_brks <- d[, fbrks(GDD)]
tm_brks <- d[, fbrks(tm)]
pp_brks <- c(seq(0, 4000, l = 50), 8200) 
ai_brks <- d[, fbrks(ai)]


# create bins ---------------------
dbins <- paste0(divs, "_bin")
fbins <- function(x, brks) findInterval(x, brks, rightmost.closed = T)
d[, (dbins):= mapply(fbins, .SD, dbrks, SIMPLIFY = F), .SDcols = divs]
d[, GDD_bin:= fbins(GDD, gdd_brks)]
d[, tm_bin:= fbins(tm, tm_brks)]
d[, pp_bin:= fbins(pp, pp_brks)]
d[, aip_bin:= fbins(ai, ai_brks)]


# load functions -------------
source("G:/My Drive/globcropdiv/Functions/my_plot.R")

## Att D Avg Temp ----------
max.area = 12000000
ibrks <- seq(0, max.area, l = 513)


fig.file = "G:/My Drive/globcropdiv/Plots/D_att_fx_tm.png"
# Delete file if it exists
if(file.exists(fig.file)) file.remove(fig.file)
png(filename = fig.file, units = 'in', width = 4.2, 
    height = 2, type = "cairo", res = 300, pointsize = 11)

par(mgp = c(0,0.9,0), oma = c(3,3,2,0), mar = c(.8,0.2,.2,0.2))
layout(t(1:3), width = c(5,5,1), height = c(1,1))
xr <- c(-15, 35)
yr <- c(0,50)

plot_frame(xr, yr, ylabel = bquote(italic(aD)))
ad_plot(d, "tm", "att_D_eco", tm_brks, dbrks$att_D_eco, add = T, ibrks = ibrks, showmax = T)
mtext("a", adj = 0.95, line = -1.8, font = 2, cex = 1.1)

plot_frame(xr, yr, yl = F)
ad_plot(d, "tm", "att_D_sdm", tm_brks, dbrks$att_D_sdm, add = T, ibrks = ibrks, showmax = T)

mtext("Temperature (°C)", side = 1, outer = T, line = 1.5, adj = 10/11 * .5)
mtext("b", adj = 0.95, line = -1.8, font = 2, cex = 1.1)

legend_image <- as.raster(matrix(viridis::mako(513), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '') 
text(x = 0.7, y = 0.95, 'Mha', cex = 1.1, adj = 0.5)
rasterImage(legend_image, 0, 0.1, 1, .85)
text(x = 1.2, y = seq(0.12,0.83, l = 5), adj = 0, 
     labels = seq(0, max.area, l = 5)/1e6, cex = 0.9)

dev.off()




















































# OLD R ----------------------------------------------------------------------

# plot functions -------------
# x = "GDD"; y = "att_D_eco"; xbrks = gdd_brks; ybrks = dbrks$att_D_eco; add = T; add_loess = T; ibrks = NULL; contr = F; minz = NULL; showmax = T

f <- function(d, x, y, xbrks, ybrks, add = F, add_loess = T, 
              ibrks = NULL, contr = F, minz = NULL, showmax = F, ncolors = 512){
  x_bin <- paste0(x, "_bin")
  y_bin <- paste0(y, "_bin")
  dd <- d[, sum(tcl), by = c(x_bin, y_bin)]
  dd <- dd[complete.cases(dd)]
  dd <- as.matrix(dd)
  m <- matrix(nrow = length(xbrks) - 1, ncol = length(ybrks) - 1)
  m[dd[, 1:2]] <- dd[,3]
  
  if(!is.null(minz)){
    m[m < minz] <- NA
  }
  
  xmp <- (xbrks[-1] - xbrks[-length(xbrks)])/2 + xbrks[-1]
  ymp <- (ybrks[-1] - ybrks[-length(ybrks)])/2 + ybrks[-1]
  
  if(is.null(ibrks)){
    ibrks <- seq(min(m, na.rm = T), max(m, na.rm = T), l = ncolors + 1)
  } 
  
  image(xmp, ymp, z = m, col = viridis::mako(ncolors, direction = -1), 
        xlab = x, ylab = y, add = add, breaks = ibrks)
  
  if(add_loess){
    dd <- as.data.frame(dd)
    names(dd) <- c("x", "y", "w")
    setDT(dd)
    dx <- data.table(x = 1:length(xmp), xmp = xmp)
    dy <- data.table(y = 1:length(ymp), ymp = ymp)
    dd <- dx[dd, on = "x"]
    dd <- dy[dd, on = "y"]
    setorder(dd, xmp, ymp)
    lo <- loess(ymp ~ xmp, dd, dd$w)
    lines(dd$xmp, predict(lo), col = "red")
  }
  
  if(contr){
    contour(xmp, ymp, m, add = T, drawlabels = F)
  }
  if(showmax) message(paste("max z value is", max(m, na.rm = T)))
}


fp <- function(xr, yr, xlabel = NULL, ylabel = NULL, xl = T, yl = T){
  plot(1, axes = F, type = "n", xlim = xr, ylim = yr, ylab = "", xlab = "")
  xn <- pmin(pretty(xr), xr[2])
  yn <- pretty(yr)
  axis(side = 1, at = xn, pos = yr[1], lwd.ticks = 0, tcl = -0.5, labels = F)
  axis(side = 1, at = xn[-length(xn)], pos = yr[1], tcl = -0.5, labels = xl)
  axis(side = 2, at = yn, pos = xr[1], tcl = -0.5, labels = yl)
  axis(side = 3, at = xn, tick = T, lwd.ticks = 0, labels = F, pos = yr[2])
  axis(side = 4, at = yn, tick = T, lwd.ticks = 0, labels = F, pos = xr[2])
  if(!is.null(xlabel)) mtext(side = 1, xlabel, line = 2)
  if(!is.null(ylabel)) mtext(side = 2, ylabel, line = 1.8, las = 0)
}

# Create plots --------------------
## GDD ----------
max.area = 11000000
ibrks <- seq(0, max.area, l = 513)


fig.file = "G:/My Drive/globcropdiv/Plots/D_fx_GDD.png"
# Delete file if it exists
if(file.exists(fig.file)) file.remove(fig.file)
png(filename = fig.file, units = 'in', width = 4.4, 
    height = 4, type = "cairo", res = 300, pointsize = 11)

par(mgp = c(0,0.9,0), oma = c(3,3,2,0), mar = c(.8,0.2,.2,0.2))
layout(matrix(c(1,2,5,3,4,5), nrow = 2, byrow = T), 
       width = c(5,5,1), height = c(1,1))
xr <- c(0,35)
yr <- c(0, 200)

fp(xr, yr, ylabel = bquote(Potential ~ italic(D)), xl = F)
f(d, "GDD", "pot_D_eco", gdd_brks, dbrks$pot_D_eco, add = T, ibrks = ibrks)

fp(xr, yr, yl = F, xl = F)
f(d, "GDD", "pot_D_sdm", gdd_brks, dbrks$pot_D_sdm, add = T, ibrks = ibrks)

yr <- c(0,50)
fp(xr, yr, ylabel = bquote(Attainable ~ italic(D)))
f(d, "GDD", "att_D_eco", gdd_brks, dbrks$att_D_eco, add = T, ibrks = ibrks)

fp(xr, yr, yl = F)
f(d, "GDD", "att_D_sdm", gdd_brks, dbrks$att_D_sdm, add = T, ibrks = ibrks)

mtext("Temperature", side = 1, outer = T, line = 1, adj = 10/11 * .5)

legend_image <- as.raster(matrix(viridis::mako(513), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '') 
text(x = 0.5, y = 0.75, 'Mha', cex = .9,)
rasterImage(legend_image, 0, 0.35, 1, .7)
text(x = 1.08, y = seq(0.35,0.7, l = 5), adj = 0, 
     labels = seq(0, max.area, l = 5)/1e6, cex = 0.8)

mtext("Ecocrop", outer = T, adj = .2)
mtext("SDM", outer = T, adj = .7)

dev.off()

## Avg Temp ----------
max.area = 13000000
ibrks <- seq(0, max.area, l = 513)


fig.file = "G:/My Drive/globcropdiv/Plots/D_fx_tm.png"
# Delete file if it exists
if(file.exists(fig.file)) file.remove(fig.file)
png(filename = fig.file, units = 'in', width = 4.4, 
    height = 4, type = "cairo", res = 300, pointsize = 11)

par(mgp = c(0,0.9,0), oma = c(3,3,2,0), mar = c(.8,0.2,.2,0.2))
layout(matrix(c(1,2,5,3,4,5), nrow = 2, byrow = T), 
       width = c(5,5,1), height = c(1,1))
xr <- c(-15, 35)
yr <- c(0, 200)

fp(xr, yr, ylabel = bquote(Potential ~ italic(D)), xl = F)
f(d, "tm", "pot_D_eco", tm_brks, dbrks$pot_D_eco, add = T, ibrks = ibrks, showmax = T)

fp(xr, yr, yl = F, xl = F)
f(d, "tm", "pot_D_sdm", tm_brks, dbrks$pot_D_sdm, add = T, ibrks = ibrks, showmax = T)

yr <- c(0,50)
fp(xr, yr, ylabel = bquote(Attainable ~ italic(D)))
f(d, "tm", "att_D_eco", tm_brks, dbrks$att_D_eco, add = T, ibrks = ibrks, showmax = T)

fp(xr, yr, yl = F)
f(d, "tm", "att_D_sdm", tm_brks, dbrks$att_D_sdm, add = T, ibrks = ibrks, showmax = T)

mtext("Temperature", side = 1, outer = T, line = 1, adj = 10/11 * .5)

legend_image <- as.raster(matrix(viridis::mako(513), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '') 
text(x = 0.5, y = 0.75, 'Mha', cex = .9,)
rasterImage(legend_image, 0, 0.35, 1, .7)
text(x = 1.08, y = seq(0.35,0.7, l = 5), adj = 0, 
     labels = seq(0, max.area, l = 5)/1e6, cex = 0.8)

mtext("Ecocrop", outer = T, adj = .2)
mtext("SDM", outer = T, adj = .7)

dev.off()


