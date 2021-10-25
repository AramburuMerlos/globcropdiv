library(terra)
library(data.table)
library(magrittr)

if(system('hostname', TRUE) == "ESP-RH-9891"){
  setwd("D:/globcropdiv/")
} else if(system('hostname', TRUE) == "LAPTOP-ST129J47"){ 
  setwd("G:/My Drive/globcropdiv/")
} # else if { ... 


fd <- function(x) exp(-sum(x * log(x), na.rm = T))

# Total Area ######################################
totcl <- rast("InData/TotalCropland.tif")

cells <- values(totcl) %>% {which(!is.na(.))}

dd <- data.table(lat = yFromCell(totcl, cells), area = extract(totcl, cells)[,1])
dd[, lat:= round(lat)]
dd <- dd[, .(area = sum(area)), by = "lat"]


### Maxent ############################
#### as is ##############################
sdm_suit <- "OutData/SDM/*_ME.tif" %>%  Sys.glob() %>%  rast()
crops <- names(sdm_suit)

# extract values, long 
l <- vector(mode = "list", length = length(crops))

for(i in 1:length(crops)){
  dt <- data.table(crop = crops[i],
                   cell = cells,
                   tcl = extract(totcl, cells)[,1],
                   suit = extract(sdm_suit[[i]], cells)[,1], 
                   lat = round(yFromCell(totcl, cells)))
  dt <- dt[!is.na(suit),]
  l[[i]] <- dt
  rm(dt)
}

d <- rbindlist(l)
rm(l); gc(reset = T)

p0_me <- d[suit == 0, .N]/length(cells)


# normalize suitabilities
d[, suit:= suit/sum(suit, na.rm = T), by = .(cell)]


# compute areas
d[, area:= suit * tcl]

# total crop areas and proportions by latitude 
d <- d[, .(area = sum(area)), by = .(lat, crop)]
d[, prop:= area/sum(area), by = .(lat)]
# potential diversity
d <- d[, .(pD_me = fd(prop)), by = .(lat)]
# join data sets
dd <- d[dd, on = "lat"]
# clean
rm(d, sdm_suit); gc(reset = T)


#### with 0 ##############################
sdm_suit <- "OutData/SDM/*_ME.tif" %>%  Sys.glob() %>%  rast()
crops <- names(sdm_suit)

# extract values, long 
l <- vector(mode = "list", length = length(crops))

for(i in 1:length(crops)){
  dt <- data.table(crop = crops[i],
                   cell = cells,
                   tcl = extract(totcl, cells)[,1],
                   suit = extract(sdm_suit[[i]], cells)[,1], 
                   lat = round(yFromCell(totcl, cells)))
  dt <- dt[!is.na(suit),]
  l[[i]] <- dt
  rm(dt)
}

d <- rbindlist(l)
rm(l); gc(reset = T)

# set suit < 1e-3 as 0
d[suit < 1e-3, suit:=0]

p0_me0 <- d[suit == 0, .N]/length(cells)

# normalize suitabilities
d[, suit:= suit/sum(suit, na.rm = T), by = .(cell)]

# compute areas
d[, area:= suit * tcl]

# total crop areas and proportions by latitude 
d <- d[, .(area = sum(area)), by = .(lat, crop)]
d[, prop:= area/sum(area), by = .(lat)]
# potential diversity
d <- d[, .(pD_me0 = fd(prop)), by = .(lat)]

# join data sets
dd <- d[dd, on = "lat"]
# clean
rm(d, sdm_suit); gc(reset = T)


### Random Forest ############################
sdm_suit <- "OutData/SDM/*_RF.tif" %>%  Sys.glob() %>%  rast()
crops <- names(sdm_suit)

# extract values, long 
l <- vector(mode = "list", length = length(crops))

for(i in 1:length(crops)){
  dt <- data.table(crop = crops[i],
                   cell = cells,
                   tcl = extract(totcl, cells)[,1],
                   suit = extract(sdm_suit[[i]], cells)[,1], 
                   lat = round(yFromCell(totcl, cells)))
  dt <- dt[!is.na(suit),]
  l[[i]] <- dt
  rm(dt)
}

d <- rbindlist(l)
rm(l); gc(reset = T)


p0_rf <- d[suit == 0, .N]/length(cells)

# normalize suitabilities
d[, suit:= suit/sum(suit, na.rm = T), by = .(cell)]

# compute areas
d[, area:= suit * tcl]

# total crop areas and proportions by latitude 
d <- d[, .(area = sum(area)), by = .(lat, crop)]
d[, prop:= area/sum(area), by = .(lat)]
# potential diversity
d <- d[, .(pD_rf = fd(prop)), by = .(lat)]
# join data sets
dd <- d[dd, on = "lat"]
# clean
rm(d, sdm_suit); gc(reset = T)




### Boosted Regression Trees ############################
sdm_suit <- "OutData/SDM/*_BRT.tif" %>%  Sys.glob() %>%  rast()
crops <- names(sdm_suit)

# extract values, long 
l <- vector(mode = "list", length = length(crops))

for(i in 1:length(crops)){
  dt <- data.table(crop = crops[i],
                   cell = cells,
                   tcl = extract(totcl, cells)[,1],
                   suit = extract(sdm_suit[[i]], cells)[,1], 
                   lat = round(yFromCell(totcl, cells)))
  dt <- dt[!is.na(suit),]
  l[[i]] <- dt
  rm(dt)
}

d <- rbindlist(l)
rm(l); gc(reset = T)

p0_brt <- d[suit == 0, .N]/length(cells)

# normalize suitabilities
d[, suit:= suit/sum(suit, na.rm = T), by = .(cell)]

# compute areas
d[, area:= suit * tcl]

# total crop areas and proportions by latitude 
d <- d[, .(area = sum(area)), by = .(lat, crop)]
d[, prop:= area/sum(area), by = .(lat)]
# potential diversity
d <- d[, .(pD_brt = fd(prop)), by = .(lat)]
# join data sets
dd <- d[dd, on = "lat"]
# clean
rm(d, sdm_suit); gc(reset = T)


dd

dd <- dd[area > 1e4, ]

## All Gamma ##################
sapply(dd, range)

lo_pD_me <- loess(pD_me ~ lat, data = dd, span = 0.2)
lo_pD_me0 <- loess(pD_me0 ~ lat, data = dd, span = 0.2)
lo_pD_rf <- loess(pD_rf ~ lat, data = dd, span = 0.2)
lo_pD_brt <- loess(pD_brt ~ lat, data = dd, span = 0.2)


png("G:/My Drive/globcropdiv/Plots/pD_fx_lat_bySDM.png", 
    height = 4, width = 6,units = "in", res = 300)

par(mfrow = c(1,1), mar = c(2,1,1,1), oma = c(3,3,0,0), xpd = NA, las = 1,
    mgp = c(2.5,.8,0), cex.lab = 1.3)

pcols <- viridis::viridis(3, alpha = 0.6)
lcols <- viridis::viridis(3)

plot(dd$lat ~ dd$pD_me0, 
     pch = 21, bg = pcols[1], cex = 0.5, cex.axis = 1.1,
     xlab = bquote(italic("pD")), 
     ylab = "latitude",
     xlim = c(1,125),
     ylim = c(-70,70))

lines(dd$lat ~ predict(lo_pD_me0),  col = lcols[1], lwd = 3)


points(dd$lat ~ dd$pD_rf, 
       pch = 21, bg = pcols[2], cex = 0.5)
lines(dd$lat ~ predict(lo_pD_rf), col = lcols[2], lwd = 3)

points(dd$lat ~ dd$pD_brt, 
       pch = 21, bg = pcols[3], cex = 0.5)
lines(dd$lat ~ predict(lo_pD_brt), col = lcols[3], lwd = 3)

legend("bottom", 
       legend = c("RF", "BRT", "ME"),
       x.intersp = 0.5, pch = 21, pt.bg = pcols[c(2,3,1)], box.lty = 0, ncol = 3,
       cex = 1.2, bty = 'n')

clip(-1, 127, -90, 90)
abline(h = c(23.4365, -23.4365), lty = 3)
abline(h = 0)

dev.off()

cor(dd$pD_me0, dd$pD_me)
