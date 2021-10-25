# 3d frequencies plot functions -------------

# x = "GDD"; y = "att_D_eco"; xbrks = gdd_brks; ybrks = dbrks$att_D_eco; add = T; add_loess = T; ibrks = NULL; contr = F; minz = NULL; showmax = T

# area density plot
ad_plot <- function(d, x, y, xbrks, ybrks, add = F, add_loess = T, 
              ibrks = NULL, contr = F, minz = NULL, showmax = F, 
              low_white = F, ncolors = 512){
  x_bin <- paste0(x, "_bin")
  y_bin <- paste0(y, "_bin")
  dd <- d[, sum(tcl), by = c(x_bin, y_bin)]
  dd <- dd[complete.cases(dd)]
  dd <- as.matrix(dd)
  m <- matrix(nrow = length(xbrks) - 1, ncol = length(ybrks) - 1)
  m[dd[, 1:2]] <- dd[,3]
  
  if(!is.null(minz)){
    m[m < minz] <- NA
  }
  
  xmp <- (xbrks[-1] - xbrks[-length(xbrks)])/2 + xbrks[-length(xbrks)]
  ymp <- (ybrks[-1] - ybrks[-length(ybrks)])/2 + ybrks[-length(ybrks)]
  
  if(is.null(ibrks)){
    ibrks <- seq(min(m, na.rm = T), max(m, na.rm = T), l = ncolors + 1)
  } 
  
  cols <- viridis::mako(ncolors, direction = -1)
  
  if(low_white){
    cols <- viridis::mako(ncolors - 1, direction = -1)
    cols <- c("#FFFFFFFF", cols)
  } else {
    cols <- viridis::mako(ncolors, direction = -1)
  }
  
  image(xmp, ymp, z = m, col = cols, 
        xlab = x, ylab = y, add = add, breaks = ibrks)
  
  if(add_loess){
    dd <- as.data.frame(dd)
    names(dd) <- c("x", "y", "w")
    setDT(dd)
    dx <- data.table(x = 1:length(xmp), xmp = xmp)
    dy <- data.table(y = 1:length(ymp), ymp = ymp)
    dd <- dx[dd, on = "x"]
    dd <- dy[dd, on = "y"]
    setorder(dd, xmp, ymp)
    lo <- loess(ymp ~ xmp, dd, dd$w)
    lines(dd$xmp, predict(lo), col = "red")
  }
  
  if(contr){
    contour(xmp, ymp, m, add = T, drawlabels = F)
  }
  if(showmax) message(paste("max z value is", max(m, na.rm = T)))
}


plot_frame <- function(xr, yr, xlabel = NULL, ylabel = NULL, xl = T, yl = T){
  plot(1, axes = F, type = "n", xlim = xr, ylim = yr, ylab = "", xlab = "")
  if(max(pretty(xr) > xr[2])){
    xn <- pmin(pretty(xr), xr[2])  
    axis(side = 1, at = xn, pos = yr[1], lwd.ticks = 0, tcl = -0.5, labels = F)
    axis(side = 1, at = xn[-length(xn)], pos = yr[1], tcl = -0.5, labels = xl)
  } else {
    xn <- pretty(xr)
    axis(side = 1, at = xn, pos = yr[1], tcl = -0.5, labels = xl)
  }

  if(max(pretty(yr) > yr[2])){
    yn <- pmin(pretty(yr), yr[2])  
    axis(side = 2, at = yn, pos = xr[1], lwd.ticks = 0, tcl = -0.5, labels = F)
    axis(side = 2, at = yn[-length(yn)], pos = xr[1], tcl = -0.5, labels = yl)
  } else {
    yn <- pretty(yr)
    axis(side = 2, at = yn, pos = xr[1], tcl = -0.5, labels = yl)
  }

  axis(side = 3, at = xn, tick = T, lwd.ticks = 0, labels = F, pos = yr[2])
  axis(side = 4, at = yn, tick = T, lwd.ticks = 0, labels = F, pos = xr[2])
  if(!is.null(xlabel)) mtext(side = 1, xlabel, line = 2)
  if(!is.null(ylabel)) mtext(side = 2, ylabel, line = 1.8, las = 0)
}
