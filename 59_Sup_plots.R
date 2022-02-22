library(data.table)
library(terra)

if(system('hostname', TRUE) %in% c("ESP-RH-9891", "LAPTOP-ST129J47")){
  setwd("D:/globcropdiv/")
} # else if { ... 


rpD <- rast("OutData/pot_D_sdm.tif")
apD <- rast("OutData/pot_D_eco.tif")

suit <- rast(Sys.glob("OutData/Ecocrop/*.tif"))
totsuit <- app(suit, sum, na.rm = T)
totsuit <- project(totsuit, rpD, mask = TRUE)


d <- data.table(cell = which(!is.na(values(totsuit))))
d[, tot.suit:= extract(totsuit, cell)]


d[, ap:= extract(apD, cell)]
d[, rp:= extract(rpD, cell)]
d[, rp_ap:= rp - ap]


set.seed(0)
ds <- d[sample(.N, 1e4)]
setorder(ds, tot.suit)
lo <- loess(rp_ap ~ tot.suit, data = ds)

fig.file = "G:/My Drive/globcropdiv/Plots/Dp_dif_fx_sum_asuit.png"
# Delete file if it exists
if(file.exists(fig.file)) file.remove(fig.file)
png(filename = fig.file, units = 'in', width = 4.5, 
    height = 4, type = "cairo", res = 300, pointsize = 11)

par(mar = c(4,4,1,1), mgp = c(3,0.8,0), las = 1)

plot(ds$tot.suit, ds$rp_ap, col = "#00000010", 
     xlab = "", ylab = "", cex.axis = 1.1)
lines(predict(lo) ~ ds$tot.suit, col = "red", lwd = 2)
mtext(bquote(Sigma ~ "absolute suitability"), side = 1, line = 2.5, cex = 1.5)
mtext(bquote(Delta[italic(pD)]), 
      side = 2, outer = T, cex = 2, las = 0, line = -1.8)

abline(h = 0)

dev.off()




