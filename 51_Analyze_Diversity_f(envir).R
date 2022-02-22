library(data.table)
library(terra)


if(system('hostname', TRUE) == "ESP-RH-9891"){
  setwd("D:/globcropdiv/")
} else if(system('hostname', TRUE) == "LAPTOP-ST129J47"){ 
  setwd("G:/My Drive/globcropdiv/")
} # else if { ... 

# upload data --------------
cl <- rast("OutData/projected/CroplandProp.tif")
cl_mask <- cl < 0.005
totcl <- cl * prod(res(cl))/10000
totcl <- mask(totcl, cl_mask, maskvalue = 1)

rGDD <- rast("InData/WorldClim/2.1/wc5min/extra/GDD.tif")
rtm <- rast("InData/WorldClim/2.1/wc5min/bioc/wc2.1_5m_bio_1.tif")
rpp <- rast("InData/WorldClim/2.1/wc5min/bioc/wc2.1_5m_bio_12.tif")
rai <- rast("InData/WorldClim/2.1/wc5min/extra/AI.tif")

rGDD <- project(rGDD, cl, mask = TRUE)
rtm <- project(rtm, cl, mask = TRUE)
rpo <- project(rpp, cl, mask = TRUE)
rai <- project(rai, cl, mask = TRUE)


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


