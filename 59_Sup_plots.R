library(data.table)
library(terra)

setwd("D:/globcropdiv/")

# one plot ######

d <- fread("OutData/CropProp_lat8.csv")
cats <- read.csv("G:/My Drive/globcropdiv/AuxData/CropCategories2.csv")
#cats <- rbind(cats, data.frame(crop="none", crop.category="None", category=12))
d <- merge(d, cats, by = "crop")
setorderv(d, "prop", order = -1L)
d <- d[prop > 0.01]
d[, x:= rep(1:5, 4)]
d[, y:= rep(4:1, each = 5)]

colors_cat <- c(brewer.pal(9, "Set1"), brewer.pal(7, "Set2")[c(1,7)])

fig.file = "G:/My Drive/globcropdiv/Plots/crop_prop_lat8.png"
# Delete file if it exists
if(file.exists(fig.file)) file.remove(fig.file)
png(filename = fig.file, units = 'in', width = 4.5, 
    height = 4, type = "cairo", res = 300, pointsize = 11)

par(mar = c(0,0,0,11), xpd = TRUE)
d[, plot(x,y, pch = 19, cex = prop * 100, col = colors_cat[category], 
         axes = F, xlab = "", ylab = "", 
         xlim = c(0, max(x)+1), ylim = c(0, max(y)+1))]

legend(legend = c(unique(cats$crop.category)),
       cex = 1, x = "right", fill = colors_cat, inset = -.7,
       title = "Crop", title.adj = 0.5)


dev.off()



# another plot ###########


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







# OLD R #############################


reco <- rast("OutData/att_D_eco.tif")
rsdm <- rast("OutData/att_D_sdm.tif")
d[, aD_eco:= extract(reco, cell)]
d[, aD_sdm:= extract(rsdm, cell)]

cor(d$aD_eco, d$aD_sdm)

ds <- d[sample(.N, 1e4),]
plot(ds$aD_eco, ds$aD_sdm, col = "#00000020", 
     xlab = "", ylab = "", cex.axis = 0.8)
abline(a = 0, b = 1)
