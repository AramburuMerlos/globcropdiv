library(terra)
library(data.table)
library(magrittr)

if(system('hostname', TRUE) == "ESP-RH-9891"){
  setwd("D:/globcropdiv/")
} else if(system('hostname', TRUE) == "LAPTOP-ST129J47"){ 
  setwd("G:/My Drive/globcropdiv/")
} # else if { ... 

# Total diversity (Dgamma) ######################################

## Actual diversity ------
totcl <- rast("InData/TotalCropland.tif")
actual <- "InData/CropAbundance/*.tif" %>%  Sys.glob() %>%  rast()
crops <- names(actual)

# extract values, long 
cells <- values(totcl) %>% {which(!is.na(.))}

dl <- vector(mode = "list", length = length(crops))

for(i in 1:length(crops)){
  dt <- data.table(crop = crops[i], 
                   cell = cells,
                   area = extract(actual[[i]], cells)[,1])
  dt <- dt[!is.na(area),]
  dt <- dt[area > 0,]
  dt[, lat:= yFromCell(actual[[i]], cell)]
  dl[[i]] <- dt
  rm(dt)
}

d <- rbindlist(dl)
rm(dl)

d[, lat:= round(lat)]

d <- d[, .(area = sum(area)), by = .(lat, crop)]

d[, prop:= area/sum(area), by = .(lat)]

# diversity function 
fd <- function(x) exp(-sum(x * log(x), na.rm = T))

dd <- d[, .(tot.area = sum(area), act_Dg = fd(prop)), by = .(lat)]




## Attainable diversity ---------
### Ecocrop --------------

alloc_eco <- fread("OutData/allocated_eco.csv")
d <- melt(alloc_eco, id.vars = "cell", measure.vars = patterns("^a."),
          variable.name = "crop", value.name = "area")
d <- d[area > 0, ]
d[, lat:= round(yFromCell(totcl, cell))]

d <- d[, .(area = sum(area)), by = .(lat, crop)]
d[, prop := area/sum(area), by = .(lat)]

d <- d[, .(att_Dg_eco = fd(prop)), by = .(lat)]
dd <- d[dd, on = "lat"]


rm(alloc_eco); gc(reset = T)



### SDM --------------------
alloc_sdm <- fread("OutData/allocated_sdm.csv")
d <- melt(alloc_sdm, id.vars = "cell", measure.vars = patterns("^a."),
          variable.name = "crop", value.name = "area")
d <- d[area > 0, ]
d[, lat:= round(yFromCell(totcl, cell))]

d <- d[, .(area = sum(area)), by = .(lat, crop)]
d[, prop := area/sum(area), by = .(lat)]

d <- d[, .(att_Dg_sdm = fd(prop)), by = .(lat)]
dd <- d[dd, on = "lat"]



rm(alloc_sdm); gc(reset = T)

## Potential Diversity ########################

### Ecocrop ############################
eco_suit <- "OutData/Ecocrop/*.tif" %>%  Sys.glob() %>%  rast()
all.equal(names(eco_suit), crops)

# extract values, long 
l <- vector(mode = "list", length = length(crops))

for(i in 1:length(crops)){
  dt <- data.table(crop = crops[i],
                   cell = cells,
                   tcl = extract(totcl, cells)[,1],
                   suit = extract(eco_suit[[i]], cells)[,1], 
                   lat = round(yFromCell(totcl, cells)))
  dt <- dt[!is.na(suit),]
  dt <- dt[suit > 0,]
  l[[i]] <- dt
  rm(dt)
}

d <- rbindlist(l)
rm(l); gc(reset = T)

# normalize suitabilities
d[, suit:= suit/sum(suit, na.rm = T), by = .(cell)]

# compute areas
d[, area:= suit * tcl]

# total crop areas by latitude 
d <- d[, .(area = sum(area)), by = .(lat, crop)]

d[, prop:= area/sum(area), by = .(lat)]

d <- d[, .(pot_Dg_eco = fd(prop)), by = .(lat)]

dd <- d[dd, on = "lat"]


### SDM ############################
sdm_suit <- "OutData/SDM/*_AVG.tif" %>%  Sys.glob() %>%  rast()
all.equal(names(sdm_suit), crops)

# extract values, long 
l <- vector(mode = "list", length = length(crops))

for(i in 1:length(crops)){
  dt <- data.table(crop = crops[i],
                   cell = cells,
                   tcl = extract(totcl, cells)[,1],
                   suit = extract(sdm_suit[[i]], cells)[,1], 
                   lat = round(yFromCell(totcl, cells)))
  dt <- dt[!is.na(suit),]
  dt <- dt[suit > 0,]
  l[[i]] <- dt
  rm(dt)
}

d <- rbindlist(l)
rm(l); gc(reset = T)

# normalize suitabilities
d[, suit:= suit/sum(suit, na.rm = T), by = .(cell)]

# compute areas
d[, area:= suit * tcl]

# total crop areas by latitude 
d <- d[, .(area = sum(area)), by = .(lat, crop)]

d[, prop:= area/sum(area), by = .(lat)]

d <- d[, .(pot_Dg_sdm = fd(prop)), by = .(lat)]

dd <- d[dd, on = "lat"]

rm(d, eco_suit, sdm_suit, actual); gc(reset = T)

setorder(dd, lat)
dd


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





## Gaps  ##################


dd[, att_Dg:= (att_Dg_eco + att_Dg_sdm)/2]
dd[, Dgap:= (1 - (act_Dg/att_Dg))*100]

lo_Dgap <- loess(Dgap ~ lat, data = dd, span = 0.2)



png("G:/My Drive/globcropdiv/Plots/Dgap_fx_lat.png", 
    height = 4, width = 4,units = "in", res = 300)

par(mfrow = c(1,1))

plot(dd$lat ~ dd$Dgap, 
     pch = 21, bg = "dark red", cex = 0.5, cex.axis = 1.1,
     xlab = bquote(italic("Dg")*"(%)"), 
     ylab = "latitude",
     xlim = c(0,80),
     ylim = c(-70,70))

lines(dd$lat ~ predict(lo_Dgap),  col = "dark red", lwd = 3)

clip(0, 80, -90, 90)
abline(h = c(23.4365, -23.4365), lty = 3)
abline(h = 0)


dd[, Mha:= tot.area/1e6]

lo_tcl <- loess(Mha ~ lat, data = dd, span = 0.2)

dev.off()


# Cropland per latitude

png("G:/My Drive/globcropdiv/Plots/tcl_fx_lat.png", 
    height = 4, width = 4,units = "in", res = 300)

plot(dd$lat ~ dd$Mha, 
     pch = 22, bg = "dark green", cex = 0.5, cex.axis = 1.1,
     xlab = "Cropland area (Mha)",
     ylab = "latitude", 
     xlim = c(0, 30),
     ylim = c(-70,70))
lines(dd$lat ~ predict(lo_tcl),  col = "dark green", lwd = 3)


clip(-2, 32, -90, 90)
abline(h = c(23.4365, -23.4365), lty = 3)
abline(h = 0)


dev.off()


