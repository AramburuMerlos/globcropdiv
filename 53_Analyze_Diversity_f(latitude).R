library(terra)
library(data.table)
library(magrittr)

if(system('hostname', TRUE) == "ESP-RH-9891"){
  setwd("D:/globcropdiv/")
} else if(system('hostname', TRUE) == "LAPTOP-ST129J47"){ 
  setwd("G:/My Drive/globcropdiv/")
} # else if { ... 

cl <- rast("OutData/projected/CroplandProp.tif")

# Total diversity (Dgamma) ######################################

# diversity function 
fd <- function(x) exp(-sum(x * log(x), na.rm = T))

# uplaod proportions 
actual <- rast("OutData/projected/ActualCropProp.tif")
pot_eco <- rast("OutData/projected/PotEcoCropProp.tif")
pot_sdm <- rast("OutData/projected/PotSDMCropProp.tif")
att_eco <- rast("OutData/projected/AttEcoCropProp.tif")
att_sdm <- rast("OutData/projected/AttSDMCropProp.tif")

crops <- names(actual)

## Actual diversity ------

# extract values, wide
d <- data.table(cell = which(!is.na(values(cl))))
d[, (crops):= extract(actual, cell)]
d[, cl_prop:= extract(cl, cell)]

# add latitude
d[, lat_proj:= yFromCell(cl, cell)]
d[, lat:= rgdal::project(cbind(0, lat_proj), crs(cl, proj=TRUE), inv=TRUE)[,2]]

# from prop to area (in ha, although units doesn't matter)
d[, totcl:= cl_prop * (prod(res(cl)))/1e4]
d[, (crops):= lapply(.SD, function(x) x*totcl), .SDcols = crops]

# reshape to long
d <- melt(d, measure.vars = crops, variable.name = "crop", value.name = "area")
d <- d[!is.na(area),]
d <- d[area > 0,]

d[, lat:= round(lat)]
d <- d[, .(area = sum(area)), by = .(lat, crop)]

d[, prop:= area/sum(area), by = .(lat)]

d8 <- d[lat == 8, ] # latitude with highest cD ~ 30
fd(d8$prop)
setorderv(d8, "area", order = -1L)
fwrite(d8, "OutData/CropProp_lat8.csv")

# compute diversity
dd <- d[, .(tot.area = sum(area), act_Dg = fd(prop)), by = .(lat)]


## Attainable diversity ------

### Ecocrop -----

# extract values, wide
d <- data.table(cell = which(!is.na(values(cl))))
d[, (crops):= extract(att_eco, cell)]
d[, cl_prop:= extract(cl, cell)]

# add latitude
d[, lat_proj:= yFromCell(cl, cell)]
d[, lat:= rgdal::project(cbind(0, lat_proj), crs(cl, proj=TRUE), inv=TRUE)[,2]]

# from prop to area (in ha, although units doesn't matter)
d[, totcl:= cl_prop * (prod(res(cl)))/1e4]
d[, (crops):= lapply(.SD, function(x) x*totcl), .SDcols = crops]

# reshape to long
d <- melt(d, measure.vars = crops, variable.name = "crop", value.name = "area")
d <- d[!is.na(area),]
d <- d[area > 0,]

d[, lat:= round(lat)]
d <- d[, .(area = sum(area)), by = .(lat, crop)]

d[, prop:= area/sum(area), by = .(lat)]

# compute diversity
d <- d[, .(tot.area = sum(area), att_Dg_eco = fd(prop)), by = .(lat)]
dd <- d[dd, on = "lat"]

gc(reset = T)


### SDM --------------------

# extract values, wide
d <- data.table(cell = which(!is.na(values(cl))))
d[, (crops):= extract(att_sdm, cell)]
d[, cl_prop:= extract(cl, cell)]

# add latitude
d[, lat_proj:= yFromCell(cl, cell)]
d[, lat:= rgdal::project(cbind(0, lat_proj), crs(cl, proj=TRUE), inv=TRUE)[,2]]

# from prop to area (in ha, although units doesn't matter)
d[, totcl:= cl_prop * (prod(res(cl)))/1e4]
d[, (crops):= lapply(.SD, function(x) x*totcl), .SDcols = crops]

# reshape to long
d <- melt(d, measure.vars = crops, variable.name = "crop", value.name = "area")
d <- d[!is.na(area),]
d <- d[area > 0,]

d[, lat:= round(lat)]
d <- d[, .(area = sum(area)), by = .(lat, crop)]

d[, prop:= area/sum(area), by = .(lat)]

# compute diversity
d <- d[, .(tot.area = sum(area), att_Dg_sdm = fd(prop)), by = .(lat)]
dd <- d[dd, on = "lat"]

gc(reset = T)


## Potential Diversity ########################

### Ecocrop -----

# extract values, wide
d <- data.table(cell = which(!is.na(values(cl))))
d[, (crops):= extract(pot_eco, cell)]
d[, cl_prop:= extract(cl, cell)]

# add latitude
d[, lat_proj:= yFromCell(cl, cell)]
d[, lat:= rgdal::project(cbind(0, lat_proj), crs(cl, proj=TRUE), inv=TRUE)[,2]]

# from prop to area (in ha, although units doesn't matter)
d[, totcl:= cl_prop * (prod(res(cl)))/1e4]
d[, (crops):= lapply(.SD, function(x) x*totcl), .SDcols = crops]

# reshape to long
d <- melt(d, measure.vars = crops, variable.name = "crop", value.name = "area")
d <- d[!is.na(area),]
d <- d[area > 0,]

d[, lat:= round(lat)]
d <- d[, .(area = sum(area)), by = .(lat, crop)]

d[, prop:= area/sum(area), by = .(lat)]

# compute diversity
d <- d[, .(tot.area = sum(area), pot_Dg_eco = fd(prop)), by = .(lat)]
dd <- d[dd, on = "lat"]

gc(reset = T)


### SDM --------------------

# extract values, wide
d <- data.table(cell = which(!is.na(values(cl))))
d[, (crops):= extract(pot_sdm, cell)]
d[, cl_prop:= extract(cl, cell)]

# add latitude
d[, lat_proj:= yFromCell(cl, cell)]
d[, lat:= rgdal::project(cbind(0, lat_proj), crs(cl, proj=TRUE), inv=TRUE)[,2]]

# from prop to area (in ha, although units doesn't matter)
d[, totcl:= cl_prop * (prod(res(cl)))/1e4]
d[, (crops):= lapply(.SD, function(x) x*totcl), .SDcols = crops]

# reshape to long
d <- melt(d, measure.vars = crops, variable.name = "crop", value.name = "area")
d <- d[!is.na(area),]
d <- d[area > 0,]

d[, lat:= round(lat)]
d <- d[, .(area = sum(area)), by = .(lat, crop)]

d[, prop:= area/sum(area), by = .(lat)]

# compute diversity
d <- d[, .(tot.area = sum(area), pot_Dg_sdm = fd(prop)), by = .(lat)]
dd <- d[dd, on = "lat"]

gc(reset = T)

setorder(dd, lat)


fwrite(dd, "OutData/DfxLat.csv")

## PLOT ##########
# Read data
dd <- fread("OutData/DfxLat.csv")
dd <- dd[tot.area > 1e5,]


## loess models ###############
lo_act_Dg <- loess(act_Dg ~ lat, data = dd, span = 0.2)
lo_att_Dg_eco <- loess(att_Dg_eco ~ lat, data = dd, span = 0.2)
lo_att_Dg_sdm <- loess(att_Dg_sdm ~ lat, data = dd, span = 0.2)
lo_pot_Dg_eco <- loess(pot_Dg_eco ~ lat, data = dd, span = 0.2)
lo_pot_Dg_sdm <- loess(pot_Dg_sdm ~ lat, data = dd, span = 0.2)



## All Gamma ##################
png("G:/My Drive/globcropdiv/Plots/Gamma_Div_fx_lat_all.png", 
    height = 4, width = 8,units = "in", res = 300)

par(mfrow = c(1,2), mar = c(2,1,1,1), oma = c(3,3,0,0), xpd = NA, las = 1,
    mgp = c(2.5,.8,0), cex.lab = 1.3)

pcols <- viridis::viridis(3, alpha = 0.6)
lcols <- viridis::viridis(3)

plot(dd$lat ~ dd$act_Dg, 
     pch = 21, bg = pcols[1], cex = 0.5, cex.axis = 1.1,
     xlab = bquote(italic("D")), 
     ylab = "latitude",
     xlim = c(1,40),
     ylim = c(-70,70))

lines(dd$lat ~ predict(lo_act_Dg),  col = lcols[1], lwd = 3)

points(dd$lat ~ dd$att_Dg_sdm, 
       pch = 21, bg = pcols[2], cex = 0.5)
lines(dd$lat ~ predict(lo_att_Dg_sdm), col = lcols[2], lwd = 3)

points(dd$lat ~ dd$att_Dg_eco, 
       pch = 21, bg = pcols[3], cex = 0.5)
lines(dd$lat ~ predict(lo_att_Dg_eco), col = lcols[3], lwd = 3)

legend("bottom", 
       legend = c(expression(italic(cD)),
                  expression("rs-"*italic(aD)),
                  expression("as-"*italic(aD))),
       x.intersp = 0.5, pch = 21, pt.bg = pcols, box.lty = 0, ncol = 3,
       cex = 1.2, bty = 'n')

mtext("a", line = -2, adj = 0.95, font = 2, cex = 1.8)

clip(-0.3, 41.3, -90, 90)
abline(h = c(23.4365, -23.4365), lty = 3)
abline(h = 0)



plot(dd$lat ~ dd$pot_Dg_sdm, 
     pch = 22, bg = pcols[2], cex = 0.5, cex.axis = 1.1,
     xlab = bquote(italic("D")),
     ylab = "", yaxt = "n",
     xlim = c(45,160),
     ylim = c(-70,70))

points(dd$lat ~ dd$pot_Dg_eco, 
       pch = 22, bg = pcols[3], cex = 0.5)

legend("bottom", ncol = 3,  col = c("white", 'black', 'black'),
       legend = c("",
                  expression("rs-"*italic(pD)),
                  expression("as-"*italic(pD))),
       x.intersp = 0.5, pch = 22, pt.bg = c('white', pcols[2:3]), box.lty = 0, 
       cex = 1.2, bty = 'n')

mtext("b", line = -2, adj = 0.95, font = 2, cex = 1.8)


clip(41, 164, -90, 90)
abline(h = c(23.4365, -23.4365), lty = 3)
abline(h = 0)
lines(dd$lat ~ predict(lo_pot_Dg_sdm),  col = lcols[2], lwd = 3)
lines(dd$lat ~ predict(lo_pot_Dg_eco), col = lcols[3], lwd = 3)



#legend("bottomright", 
#       legend = c("Potential (SDM)", "Potential (Ecocrop)"), 
#       box.lty = 0, cex = 0.7, pch = rep(NA, 3), lty = 3, col = lcols[2:3])
#
dev.off()

