library(magrittr)
library(data.table)
library(terra)

if(!grepl("globcropdiv$", getwd())){
  if(system('hostname', TRUE) %in% c("ESP-RH-9891", "LAPTOP-ST129J47")){ 
    setwd("G:/My Drive/globcropdiv/")
  } # else if { ... 
}


# Da by crop ##############
# import rasters
totcl <- rast("InData/TotalCropland.tif")
Da <- rast("OutData/Da.tif")
actual <- "InData/CropAbundance/*.tif" %>%  Sys.glob() %>%  rast()
crops <- names(actual)

# extract values 
d <- values(totcl) %>% {data.table(cell = which(!is.na(.)))}
d[, (crops):= extract(actual, cell)]
d[, Da:= extract(Da, cell)]

# calculate average Da and total crop area
Da_crop <- paste0("Da_", crops)
d[, (Da_crop):= lapply(.SD, function(x) x/sum(x) * Da), .SDcols = crops]

dc <- data.table(crop = crops,
                 crop_area = colSums(d[, ..crops]), 
                 crop_Da = colSums(d[, ..Da_crop]))

cropcat <- fread("AuxData/CropCategories.csv")

dc <- merge(dc, cropcat, on = "crop")

dc[, `:=`(area = sum(crop_area), n_cat = .N), by = SPAM_Name]

ds <- dc[, .(area = mean(area/n_cat), 
             Da = sum(crop_Da * crop_area/area)), 
         by = "SPAM_Name"]

setorderv(ds, "area", order = -1L)

ds <- ds[!(SPAM_Name %like% "other" | SPAM_Name %like% "rest")]


df <- data.frame(x = log2(ds$area/sum(ds$area)), y = ds$Da)
sapply(df, range)

dir.create("Plots")

## complete plot -----
# plot
{ # run this line to save fig to folder
  fig.file = "Plots/CropDa_f(CropArea).pdf"
  #Delete file if it exists
  if(file.exists(fig.file)) file.remove(fig.file)
  pdf(file = fig.file, width = 11.4/2.54, height = (11.4/2.54 - 1.5))
  { # run here to get plot in R
    par(mfrow = c(1,1), mai = c(0.8,0.5,0,1.8), oma = c(0,0,0,0), las = 1, xpd = T)
    xr = c(-10, -2)
    yr = c(2, 12)
    xn <- length(df$x)
    cols <- viridisLite::viridis(xn)
    pt <- rep(c(21:25,0:9),10)[1:xn]
    plot(1, axes = F, type = "n", xlim = xr, ylim = yr, xlab = "", ylab = "")
    axis(side = 1, at = pretty(xr), pos = yr[1], lwd = 1, cex.axis = .85, 
         tcl = -0.4, labels = round(2^(pretty(xr))*100, 2), cex.axis = 0.7, mgp = c(0,0.7,0))
    axis(side = 2, at = pretty(yr), pos = xr[1], lwd = 1, cex.axis = 0.85, tcl = -0.4, cex.axis = 0.7)
    axis(side = 3, at = pretty(xr), tick = T, lwd.ticks = 0, labels = F, pos = yr[2], lwd = 1)
    axis(side = 4, at = pretty(yr), tick = T, lwd.ticks = 0, labels = F, pos = xr[2], lwd = 1)
    
    points(df$x, df$y, pch = pt, cex = 0.8, col = ifelse(pt < 20, cols,'black'), bg = cols)
    
    mtext(side = 1, 'Crop area (%)', line = 1.5, cex = 1)
    mtext(side = 2, 'Actual Diversity', line = 1.5, cex = 1, las = 0)
    legend(x = max(xr)+(xr[2]-xr[1])*.05, y = max(yr)*1.05, 
           legend = ds$SPAM_Name[1:24], col = ifelse(pt<20,cols,'black')[1:24], bty = 'n', 
           pch = pt[1:24], pt.bg = cols[1:24], cex = 0.58)
    legend(x = max(xr)+(xr[2]-xr[1])*.48,  y = max(yr)*1.05, 
           legend = ds$SPAM_Name[25:xn], col = ifelse(pt<20,cols,'black')[25:xn], bty = 'n', 
           pch = pt[25:xn], pt.bg = cols[25:xn], cex = 0.58)
#    clip(xr[1],xr[2], yr[1], yr[2])
#    lines(fitted, col = 'red', lty = 2)
#    lines(fitted2$x,fitted2$y, col = 'blue', lty = 4)
#    lines(tmax$lp, tmax$m, col = 'grey50', lty = 3)
#    legend(x = xr[1], y = 2, cex = 0.6, col = c("red", "blue", "grey50"), lty = c(2,4,3), bty = 'n',
#           legend = c("Empirical", "Empirical (1, 100)", "Theoretical Maximum"))
    
    #text(xr[1]+(xr[2]-xr[1])*.05,yr[1]+(yr[2]-yr[1])*.15, pos = 4, label = "p-value < 0.001")
  }
  dev.off()
}

## legend ------
legend <- 'Crop categories as reported in SPAM. For agregated crop categories (i.e. tropical fruit, temperate fruit and vegetables), the crop area is the average crop area of their crop components (rather than the sum) and the diversity is the area-weighted average diversity. All "other" crop categories and "rest of crops" were ommited.' 

writeLines(legend, "Plots/CropDa_f(CropArea)_legend.txt")

# Dp by crop ##############
# import data
d <- fread("OutData/allocated_sdm.csv")

acols <- names(d)[names(d) %like% "^a\\."]
crops <-  gsub("^a.", "", acols)

# extract diversity values 
Dp <- rast("OutData/Dp_sdm.tif")
d[, Dp:= extract(Dp, cell)]

# calculate average Dp and total crop area
Dp_crop <- paste0("Dp_", crops)
d[, (Dp_crop):= lapply(.SD, function(x) x/sum(x) * Dp), .SDcols = acols]

dc <- data.table(crop = crops,
                 crop_area = colSums(d[, ..acols]), 
                 crop_Dp = colSums(d[, ..Dp_crop]))

cropcat <- fread("AuxData/CropCategories.csv")

dc <- merge(dc, cropcat, on = "crop")

dc[, `:=`(area = sum(crop_area), n_cat = .N), by = SPAM_Name]

ds <- dc[, .(area = mean(area/n_cat), 
             Dp = sum(crop_Dp * crop_area/area)), 
         by = "SPAM_Name"]

setorderv(ds, "area", order = -1L)
ds
plot(ds$Dp, log(ds$area))
ds <- ds[!(SPAM_Name %like% "other" | SPAM_Name %like% "rest")]


df <- data.frame(x = log2(ds$area/sum(ds$area)), y = ds$Dp)
sapply(df, range)

plot(df$x, df$y)


# complete plot ####
# plot
{ # run this line to save fig to folder
  fig.file = "Plots/CropDp_f(CropArea).pdf"
  #Delete file if it exists
  if(file.exists(fig.file)) file.remove(fig.file)
  pdf(file = fig.file, width = 11.4/2.54, height = (11.4/2.54 - 1.5))
  { # run here to get plot in R
    par(mfrow = c(1,1), mai = c(0.8,0.5,0,1.8), oma = c(0,0,0,0), las = 1, xpd = T)
    xr = c(-10, -2)
    yr = c(14, 28)
    xn <- length(df$x)
    cols <- viridisLite::viridis(xn)
    pt <- rep(c(21:25,0:9),10)[1:xn]
    plot(1, axes = F, type = "n", xlim = xr, ylim = yr, xlab = "", ylab = "")
    axis(side = 1, at = pretty(xr), pos = yr[1], lwd = 1, cex.axis = .85, 
         tcl = -0.4, labels = round(2^(pretty(xr))*100, 2), cex.axis = 0.7, mgp = c(0,0.7,0))
    axis(side = 2, at = pretty(yr), pos = xr[1], lwd = 1, cex.axis = 0.85, tcl = -0.4, cex.axis = 0.7)
    axis(side = 3, at = pretty(xr), tick = T, lwd.ticks = 0, labels = F, pos = yr[2], lwd = 1)
    axis(side = 4, at = pretty(yr), tick = T, lwd.ticks = 0, labels = F, pos = xr[2], lwd = 1)
    
    points(df$x, df$y, pch = pt, cex = 0.8, col = ifelse(pt < 20, cols,'black'), bg = cols)
    
    mtext(side = 1, 'Crop area (%)', line = 1.5, cex = 1)
    mtext(side = 2, 'Potential Diversity', line = 1.5, cex = 1, las = 0)
    legend(x = max(xr)+(xr[2]-xr[1])*.05, y = max(yr)*1.05, 
           legend = ds$SPAM_Name[1:24], col = ifelse(pt<20,cols,'black')[1:24], bty = 'n', 
           pch = pt[1:24], pt.bg = cols[1:24], cex = 0.58)
    legend(x = max(xr)+(xr[2]-xr[1])*.48,  y = max(yr)*1.05, 
           legend = ds$SPAM_Name[25:xn], col = ifelse(pt<20,cols,'black')[25:xn], bty = 'n', 
           pch = pt[25:xn], pt.bg = cols[25:xn], cex = 0.58)
    #    clip(xr[1],xr[2], yr[1], yr[2])
    #    lines(fitted, col = 'red', lty = 2)
    #    lines(fitted2$x,fitted2$y, col = 'blue', lty = 4)
    #    lines(tmax$lp, tmax$m, col = 'grey50', lty = 3)
    #    legend(x = xr[1], y = 2, cex = 0.6, col = c("red", "blue", "grey50"), lty = c(2,4,3), bty = 'n',
    #           legend = c("Empirical", "Empirical (1, 100)", "Theoretical Maximum"))
    
    #text(xr[1]+(xr[2]-xr[1])*.05,yr[1]+(yr[2]-yr[1])*.15, pos = 4, label = "p-value < 0.001")
  }
  dev.off()
}


