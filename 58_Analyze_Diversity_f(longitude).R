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
  dt[, lon:= xFromCell(actual[[i]], cell)]
  dl[[i]] <- dt
  rm(dt)
}

d <- rbindlist(dl)
rm(dl)

d[, lon:= round(lon)]

d <- d[, .(area = sum(area)), by = .(lon, crop)]

d[, prop:= area/sum(area), by = .(lon)]

# diversity function 
fd <- function(x) exp(-sum(x * log(x), na.rm = T))

dd <- d[, .(tot.area = sum(area), act_Dg = fd(prop)), by = .(lon)]
#dd <- dd[tot.area > 1e5,]



## Attainable diversity ---------
### Ecocrop --------------

alloc_eco <- fread("OutData/allocated_eco.csv")
d <- melt(alloc_eco, id.vars = "cell", measure.vars = patterns("^a."),
          variable.name = "crop", value.name = "area")
d <- d[area > 0, ]
d[, lon:= round(xFromCell(totcl, cell))]

d <- d[, .(area = sum(area)), by = .(lon, crop)]
d[, prop := area/sum(area), by = .(lon)]

d <- d[, .(att_Dg_eco = fd(prop)), by = .(lon)]
dd <- d[dd, on = "lon"]


rm(alloc_eco); gc(reset = T)



### SDM --------------------
alloc_sdm <- fread("OutData/allocated_sdm.csv")
d <- melt(alloc_sdm, id.vars = "cell", measure.vars = patterns("^a."),
          variable.name = "crop", value.name = "area")
d <- d[area > 0, ]
d[, lon:= round(xFromCell(totcl, cell))]

d <- d[, .(area = sum(area)), by = .(lon, crop)]
d[, prop := area/sum(area), by = .(lon)]

d <- d[, .(att_Dg_sdm = fd(prop)), by = .(lon)]
dd <- d[dd, on = "lon"]



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
                   lon = round(xFromCell(totcl, cells)))
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

# total crop areas by longitude
d <- d[, .(area = sum(area)), by = .(lon, crop)]

d[, prop:= area/sum(area), by = .(lon)]

d <- d[, .(pot_Dg_eco = fd(prop)), by = .(lon)]

dd <- d[dd, on = "lon"]


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
                   lon = round(xFromCell(totcl, cells)))
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

# total crop areas by longitude 
d <- d[, .(area = sum(area)), by = .(lon, crop)]

d[, prop:= area/sum(area), by = .(lon)]

d <- d[, .(pot_Dg_sdm = fd(prop)), by = .(lon)]

dd <- d[dd, on = "lon"]

rm(d, eco_suit, sdm_suit, actual); gc(reset = T)

setorder(dd, lon)
dd


fwrite(dd, "OutData/DfxLon.csv")




## PLOT ##########
# Read data
dd <- fread("OutData/DfxLon.csv")

dd <- dd[tot.area > 1e5]

## loess models ###############
lo_act_Dg <- loess(act_Dg ~ lon, data = dd, span = 0.2)
lo_att_Dg_eco <- loess(att_Dg_eco ~ lon, data = dd, span = 0.2)
lo_att_Dg_sdm <- loess(att_Dg_sdm ~ lon, data = dd, span = 0.2)
lo_pot_Dg_eco <- loess(pot_Dg_eco ~ lon, data = dd, span = 0.2)
lo_pot_Dg_sdm <- loess(pot_Dg_sdm ~ lon, data = dd, span = 0.2)











png("G:/My Drive/globcropdiv/Plots/Gamma_Div_fx_lon_all.png", 
    height = 8, width = 6,units = "in", res = 300)

par(mfrow = c(2,1), mar = c(3,1,2,0), oma = c(5,4,0,6), xpd = NA, las = 1,
    mgp = c(2.5,.8,0), cex.axis = 1.3)

pcols <- viridis::viridis(3, alpha = 0.6)
lcols <- viridis::viridis(3)


plot(dd$lon, dd$pot_Dg_sdm, 
     pch = 22, bg = pcols[2], cex = 0.5,
     ylab = "",
     xlab = "", xaxt = "n",
     ylim = c(45,160),
     xlim = c(-125, 155))

points(dd$lon, dd$pot_Dg_eco, 
       pch = 22, bg = pcols[3], cex = 0.5)

mtext("a", line = -2, adj = 0.95, font = 2, cex = 1.8)


lines(dd$lon, predict(lo_pot_Dg_sdm),  col = lcols[2], lwd = 3)
lines(dd$lon, predict(lo_pot_Dg_eco), col = lcols[3], lwd = 3)


legend(160, 125, col = 'black',
       legend = c(expression("rs-"*italic(pD)),
                  expression("as-"*italic(pD))),
       pch = 22, pt.bg = c(pcols[c(2,3)]), box.lty = 0, 
       cex = 1.5, bty = 'n')




plot(dd$lon, dd$act_Dg, 
     pch = 21, bg = pcols[1], cex = 0.5, 
     ylab = "", 
     xlab = "",
     xlim = c(-125, 155),
     ylim = c(1, 40))

lines(dd$lon, predict(lo_act_Dg),  col = lcols[1], lwd = 3)

points(dd$lon, dd$att_Dg_sdm, 
       pch = 21, bg = pcols[2], cex = 0.5)
lines(dd$lon, predict(lo_att_Dg_sdm), col = lcols[2], lwd = 3)

points(dd$lon, dd$att_Dg_eco, 
       pch = 21, bg = pcols[3], cex = 0.5)
lines(dd$lon, predict(lo_att_Dg_eco), col = lcols[3], lwd = 3)

mtext("b", line = -2, adj = 0.95, font = 2, cex = 1.8)

mtext("longitude", outer = T, side = 1, adj = 0.5, cex = 2, line = 0)
mtext(bquote(italic(Diversity)), outer = T, side = 2, 
      adj = 0.5, cex = 2, line = 1.7, las = 3)

legend(160, 30, col = 'black',
       legend = c(expression("rs-"*italic(aD)),
                  expression("as-"*italic(aD)),
                  expression(italic(cD))),
       pch = 21, pt.bg = c(pcols[c(2,3,1)]), box.lty = 0, 
       cex = 1.5, bty = 'n')

w <- geodata::world(path = "InData")
w <- crop(w, ext(c(-145, 175, - 60, 80)))


par(mfrow = c(1,1), oma = c(0,0,0,0), new = T)
plot(w, border = "#00000050", mar = c(10,4.4,5,6.8),
     xlim = c(-145, 155), ylim = c(-90, 90 ), axes = FALSE)


dev.off()






