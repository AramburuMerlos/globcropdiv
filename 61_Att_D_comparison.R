library(data.table)
library(terra)


if(system('hostname', TRUE) == "ESP-RH-9891"){
  setwd("D:/globcropdiv/")
} else if(system('hostname', TRUE) == "LAPTOP-ST129J47"){ 
  setwd("G:/My Drive/globcropdiv/")
} # else if { ... 

# upload data --------------

totcl <- rast("InData/TotalCropland.tif")

# extract values
d <- data.table(cells = which(!is.na(values(totcl))))
d[, tcl:= extract(totcl, cells)[,1]]

rdivs <- rast(Sys.glob("OutData/*_D*tif"))
divs <- names(rdivs)
d[, (divs):= extract(rdivs, cells)]

d
aD <- c("att_D_eco", "att_D_sdm")
pD <- gsub("att", "pot", aD)

# create breaks ---------------------
aD_brks <- 0:50
pD_brks <- seq(0,175, length.out = 51)

# create bins ---------------------
aD_bins <- paste0(aD, "_bin")
pD_bins <- paste0(pD, "_bin")

fbins <- function(x, brks) findInterval(x, brks, rightmost.closed = T)

d[, (aD_bins):= lapply(.SD, fbins, aD_brks), .SDcols = aD]
d[, (pD_bins):= lapply(.SD, fbins, pD_brks), .SDcols = pD]
d

# compute stats -------------------
fwrmse <- function(x, y, w, na.rm = F){
  w <- w/sum(w, na.rm = na.rm)
  sqrt(sum((x - y)^2 * w, na.rm = na.rm))
  
}

aD_wrmse <- d[, fwrmse(att_D_eco, att_D_sdm, tcl, na.rm = T)]
pD_wrmse <- d[, fwrmse(pot_D_eco, pot_D_sdm, tcl, na.rm = T)]

aD_r <- d[, cor(att_D_eco, att_D_sdm, use = "complete.obs")]
pD_r <- d[, cor(pot_D_eco, pot_D_sdm, use = "complete.obs")]

max.area = 16000000
ibrks <- seq(0, max.area, l = 513)

# Create plots --------------------
source("G:/My Drive/globcropdiv/Functions/my_plot.R")

fig.file = "G:/My Drive/globcropdiv/Plots/D_comp.png"
# Delete file if it exists
if(file.exists(fig.file)) file.remove(fig.file)

png(fig.file, 6.6, 3, unit = "in", res = 300)
par(mgp = c(0,.8,0), cex.axis = 1.2, mar = c(4,4,1.5,1.5))

layout(t(c(1:3)), width = c(5,5,1), height = c(1,1))

xr <- yr <- range(pD_brks)
plot_frame(xr, yr, xl = T, 
           xlabel = bquote("as-"*italic(pD)),
           ylabel = bquote("rs-"*italic(pD)))

ad_plot(d, "pot_D_eco", "pot_D_sdm", pD_brks, pD_brks, ibrks = ibrks,
        add = T, add_loess = F, low_white = T, showmax = T)
clip(0,175,0,175)
abline(a=0, b=1, col = "red")


text(100, 25, paste0("wRMSE = ", round(pD_wrmse, 2)), cex = 1.1, adj = 0)
text(100, 10, bquote(italic(r)~"="~.(round(pD_r,2))), cex = 1.1, adj = 0)
mtext("a", adj = 0.95, line = 0, font = 2, cex = 1.1)



xr <- yr <- range(aD_brks)

plot_frame(xr, yr, xl = T,
           xlabel = bquote("as-"*italic(aD)),
           ylabel = bquote("rs-"*italic(aD)))

ad_plot(d, "att_D_eco", "att_D_sdm", aD_brks, aD_brks, ibrks = ibrks,
        add = T, add_loess = F, low_white = T, showmax = T)
clip(0,50,0,50)
abline(a=0, b=1, col = "red")


text(30, 7, paste0("wRMSE = ", round(aD_wrmse, 2)), cex = 1.1, adj = 0)
text(30, 3.5, bquote(italic(r)~"="~.(round(aD_r,2))), cex = 1.1, adj = 0)
mtext("b", adj = 0.95, line = 0, font = 2, cex = 1.1)



legend_image <- as.raster(matrix(c(viridis::mako(512), "#FFFFFFFF"), ncol=1))
par(mar = c(2,0,0,0))
plot(c(0,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '') 
text(x = 0.5, y = 0.95, 'Mha', cex = 1.2,)
rasterImage(legend_image, 0, 0.1, 1, .9)
text(x = 1.2, y = seq(0.1,0.9, l = 5), adj = 0, 
     labels = seq(0, max.area, l = 5)/1e6, cex = 1.1)

dev.off()


