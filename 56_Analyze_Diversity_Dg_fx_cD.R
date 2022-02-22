library(data.table)
library(terra)
library(quantreg)

if(system('hostname', TRUE) %in% c("ESP-RH-9891", "LAPTOP-ST129J47")){
  setwd("D:/globcropdiv/")
} # else if { ... 


pcl <- rast("OutData/projected/CroplandProp.tif")
cl_mask <- pcl < 0.005
totcl <- pcl * prod(res(pcl))/1e4
totcl <- mask(totcl, cl_mask, maskvalue = 1)


cD <- rast("OutData/act_D.tif")
rs_aD <- rast("OutData/att_D_sdm.tif")
as_aD <- rast("OutData/att_D_eco.tif")

d <- data.table(cells = which(!is.na(values(totcl))))
d[, tcl:= extract(totcl, cells)[,1]]
d[, cD:= extract(cD, cells)[,1]]
d[, rs_aD:= extract(rs_aD, cells)[,1]]
d[, as_aD:= extract(as_aD, cells)[,1]]
d[is.na(as_aD), as_aD:= 1]

d[, aD:= (rs_aD + as_aD)/2]

d[, Dg:= (aD - cD)/aD * 100]

sapply(d, range)


d[Dg < 0, Dg:= 0]



# create breaks ---------------------
cD_brks <- seq(1, 50, 0.5)
Dg_brks <- 0:100

# create bins ---------------------
fbins <- function(x, brks) findInterval(x, brks, rightmost.closed = T)
d[, cD_bin:= fbins(cD, cD_brks)]
d[, Dg_bin:= fbins(Dg, Dg_brks)]


# quantil regressions ----------------
# instead of fitting with the whole data (700k rows), do it for bins
dd <- d[, .(area = sum(tcl)), by = c("cD_bin", "Dg_bin")]
dd[, w:= area/sum(area)]
dd[, cD_l:= cD_brks[cD_bin]]
dd[, cD_h:= cD_brks[cD_bin + 1]]
dd[, Dg_l:= Dg_brks[Dg_bin]]
dd[, Dg_h:= Dg_brks[Dg_bin + 1]]
dd[, cD_m:= (cD_l + cD_h)/2]
dd[, Dg_m:= (Dg_l + Dg_h)/2]

q01 <- rq(Dg_l ~ cD_l, tau = .01, data = dd, weights = w)
q25 <- rq(Dg_l ~ cD_l, tau = .25, data = dd, weights = w)
q50 <- rq(Dg_m ~ cD_m, tau = .50, data = dd, weights = w)
q75 <- rq(Dg_h ~ cD_h, tau = .75, data = dd, weights = w)
q99 <- rq(Dg_h ~ cD_h, tau = .99, data = dd, weights = w)


# distribution -----------
setorder(d, "cD")
d[, cf:= floor(cumsum(tcl)/sum(tcl) * 10)]
cD_rug <- d[, min(cD), by = "cf"][cf > 0 & cf < 10, V1]

setorder(d, "Dg")
d[, cf:= floor(cumsum(tcl)/sum(tcl) * 10)]
Dg_rug <- d[, min(Dg), by = "cf"][cf > 0 & cf < 10, V1]
d[cf == 9]


# load functions -------------
source("G:/My Drive/globcropdiv/Functions/my_plot.R")

max.area = 12000000
ibrks <- seq(0, max.area, l = 513)


fig.file = "G:/My Drive/globcropdiv/Plots/Dg_fx_cD.png"
# Delete file if it exists
if(file.exists(fig.file)) file.remove(fig.file)
png(filename = fig.file, units = 'in', width = 3.6, 
    height = 3, type = "cairo", res = 300, pointsize = 11)

xr <- c(0, 35)
yr <- c(0, 100)

par(mgp = c(0,0.9,0), oma = c(3,3,2,0), mar = c(.8,0.2,.2,0.2), las = 1)
layout(t(1:2), width = c(10,2), height = c(1,1))
plot_frame(xr, yr, 
           xlabel = bquote(italic(cD)), 
           ylabel = bquote(italic(D)*g~"(%)"))
rug(cD_rug, ticksize = 0.06, side = 3)
rug(Dg_rug, ticksize = 0.07, side = 4)





clip(1,35,0,100)
ad_plot(d, "cD", "Dg", cD_brks, Dg_brks, add = T, ibrks = ibrks, add_loess = F, 
        showmax = T)
abline(reg = q01, col = "red", lwd = .5, lty = 3)
abline(reg = q25, col = "red", lwd = .5, lty = 2)
abline(reg = q50, col = "red", lwd = .5, lty = 1)
abline(reg = q75, col = "red", lwd = .5, lty = 2)
abline(reg = q99, col = "red", lwd = .5, lty = 3)



legend_image <- as.raster(matrix(viridis::mako(513), ncol=1))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '') 
text(x = 0.7, y = 0.95, 'Mha', cex = 1.1, adj = 0.5)
rasterImage(legend_image, 0, 0.1, 1, .85)
text(x = 1.2, y = seq(0.12,0.83, l = 5), adj = 0, 
     labels = seq(0, max.area, l = 5)/1e6, cex = 0.9)

dev.off()




