library(terra)
library(data.table)
library(magrittr)

if(system('hostname', TRUE) %in% c("ESP-RH-9891", "LAPTOP-ST129J47")){
  setwd("D:/globcropdiv/")
} # else if { ... 


cl <- rast("OutData/projected/CroplandProp.tif")

# Total diversity (Dgamma) ######################################

# diversity function 
fd <- function(x) exp(-sum(x * log(x), na.rm = T))

# uplaod proportions 
actual <- rast("OutData/projected/ActualCropProp.tif")

crops <- names(actual)

## Actual diversity ------

# extract values, wide
d <- data.table(cell = which(!is.na(values(cl))))
d[, (crops):= extract(actual, cell)]
d[, cl_prop:= extract(cl, cell)]

# add lon lat
d[, lon_proj:= xFromCell(cl, cell)]
d[, lat_proj:= yFromCell(cl, cell)]
xy <- d[, rgdal::project(cbind(lon_proj, lat_proj), crs(cl, proj=TRUE), inv=TRUE)]
d[, lon:= xy[,1]]
d[, lat:= xy[,2]]

# define west
d[, west:= lon < -30]

# from prop to area (in ha, although units doesn't matter)
d[, totcl:= cl_prop * (prod(res(cl)))/1e4]
d[, (crops):= lapply(.SD, function(x) x*totcl), .SDcols = crops]

# reshape to long
d <- melt(d, measure.vars = crops, variable.name = "crop", value.name = "area")
d <- d[!is.na(area),]
d <- d[area > 0,]

d[, lat:= round(lat)]
d <- d[, .(area = sum(area)), by = .(lat, west, crop)]

d[, prop:= area/sum(area), by = .(lat, west)]

# compute diversity
dd <- d[, .(tot.area = sum(area), cD = fd(prop)), by = .(lat, west)]

rm(d); gc(reset = T)

setorder(dd, west, lat)

fwrite(dd, "OutData/DfxLat_w.csv")


## PLOT ##########
# Read data
dd <- fread("OutData/DfxLat_w.csv")

hist(log10(dd$tot.area))
dd <- dd[tot.area > 5e4,]


## loess models ###############
lo_cDw <- loess(cD ~ lat, data = dd[(west)], span = 0.3)
lo_cDe <- loess(cD ~ lat, data = dd[!(west)], span = 0.3)

## All Gamma ##################
png("G:/My Drive/globcropdiv/Plots/cD_lat_we.png", 
    height = 5, width = 5,units = "in", res = 300)

par(mfrow = c(1,1), mar = c(4,4,1,1), xpd = NA, las = 1,
    mgp = c(2.5,.8,0), cex.lab = 1.3)

pcols <- viridis::viridis(2, begin = 0.1, end = 0.8, alpha = 0.6)
lcols <- viridis::viridis(2, begin = 0.1, end = 0.8)

dd[(west), .(plot(lat ~ cD, 
                  pch = 21, bg = pcols[1], cex = 0.7, cex.axis = 1.1,
                  xlab = bquote(italic("cD")), 
                  ylab = "latitude",
                  xlim = c(0,30),
                  ylim = c(-70,70)),
             lines(lat ~ predict(lo_cDw), col = lcols[1], lwd = 3))]

dd[!(west), .(points(lat ~ cD, pch = 21, bg = pcols[2], cex = 0.7),
              lines(lat ~ predict(lo_cDe), col = lcols[2], lwd = 3))]

clip(-0.8, 31.3, -90, 90)
abline(h = c(23.4365, -23.4365), lty = 3)
abline(h = 0)



legend("bottomright", legend = c("west", "east"),
       lwd = 3, col = lcols, box.lty = 0, x.intersp = 0.5,
       cex = 1.2, bty = 'n')

legend("bottomright", legend = c("west", "east"),
       pch = 21, pt.bg = pcols, box.lty = 0, x.intersp = 1.5,
       cex = 1.2, bty = 'n')

dev.off()

