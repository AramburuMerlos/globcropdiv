library(magrittr)
library(data.table)
library(terra)

if(system('hostname', TRUE) == "ESP-RH-9891"){
  setwd("D:/globcropdiv/")
} else if(system('hostname', TRUE) == "LAPTOP-ST129J47"){ 
  setwd("G:/My Drive/globcropdiv/")
} # else if { ... 

# Data prep --------

# import actual area raster layers
totcl <- rast("InData/TotalCropland.tif")
actual <- "InData/CropAbundance/*.tif" %>%  Sys.glob() %>%  rast()
crops <- names(actual)

# extract values 
d <- values(totcl) %>% 
  {data.table(cell = which(!is.na(.)), tcl = .[!is.na(.)])}
d[, (crops):= extract(actual, cell)]

# area to fraction
for(j in crops) set(d, j = j, value = d[[j]]/d$tcl)

# import allocated area data
d_int <- fread("OutData/allocated_int.csv")

acols <- names(d_int)[names(d_int) %like% "^a\\."]
scols <- names(d_int)[names(d_int) %like% "^s\\."]


# allocated area to fraction
for(j in acols) set(d_int, j = j, value = d_int[[j]]/d_int$tcl)

# JOIN DT
setkey(d, cell)
setkey(d_int, cell)
d <- d_int[d]

rm(d_int)
gc(reset = T)
# some cells are not in d_int b/c of being estimated as unsuitable for all crops
# change NA to 0
for(j in names(d)) set(d, which(is.na(d[[j]])), j, 0)


## SPAM categories ----
cropcat <- fread("AuxData/CropCategories.csv")

# crops to be aggregated
cropcat[, nsp:= .N, by = SPAM_Code]

agcrops <- unique(cropcat$SPAM_Name[cropcat$nsp > 1])

for(j in agcrops){
  cjs = cropcat[SPAM_Name == j, crop]
  d[, (j):= Reduce(`+`, .SD), .SDcols = cjs]
  
  aj <- paste0("a.", j)
  acjs <- paste0("a.", cjs)
  d[, (aj):= Reduce(`+`, .SD), .SDcols = acjs]

  sj <- paste0("s.", j)
  scjs <- paste0("s.", cjs)
  d[, (sj):= Reduce(`+`, .SD), .SDcols = scjs]
} 

d[, (cropcat[SPAM_Name %in% agcrops, crop]):= NULL]
d[, paste0("a.", (cropcat[SPAM_Name %in% agcrops, crop])):= NULL]
d[, paste0("s.", (cropcat[SPAM_Name %in% agcrops, crop])):= NULL]

acols <- names(d)[names(d) %like% "^a\\."]
scols <- names(d)[names(d) %like% "^s\\."]
crops <- gsub("^a\\.", "", acols)

# crop area 
crop_area <- colSums(d[, lapply(.SD, function(x) x * tcl), .SDcols = crops])

# crop niche area 
crop_niche <-  colSums(d[, lapply(.SD, function(x) x * tcl), .SDcols = scols])
names(crop_niche) <- gsub("^s\\.", "", names(crop_niche))
  

# when adding crop fractions, some might turn out to be (very slightly) > 1
for(j in c(acols, crops)) d[d[[j]] > 1, (j):= 1]


# Create Bins ------
# for each crop, no bin should have an expected count of less than 5 cells. 
# Thus, crops might be different for different crops, merging bins when necessary

# list of bins
lbins <- vector(mode = "list", length = length(crops))
names(lbins) <- crops

# vector with chi square values
chisqv <- vector(mode = "numeric", length = length(crops))
names(chisqv) <- crops

nbreaks <- floor(2*d[,.N] ^(2/5))
rbreaks <- round(d[,.N]/nbreaks * 1:nbreaks)

# variable bin widths: equiprobable (expected value) bins
for(i in 1:length(crops)){
  
  oa <- crops[i]
  ea <- acols[i]

  setorderv(d, acols[i])
  d[, cs:= cumsum(.SD), .SDcols = ea]
  breaks <- c(0, d[[ea]][rbreaks])
  
  # be sure to cover all possible values so no observed value is left outside
  breaks[length(breaks)] <- 1
  # avoid repeated breaks (e.g.at 0)
  breaks <- unique(breaks)
  
  # calculate frequencies
  d[, obs_bin:= cut(d[[oa]], breaks = breaks, include.lowest = TRUE)]
  dt_obs <- d[, .N, by = obs_bin]

  d[, exp_bin:= cut(d[[ea]], breaks = breaks, include.lowest = TRUE)]
  dt_exp <- d[, .N, by = exp_bin]

  dt <- merge(dt_obs, dt_exp, by.x = "obs_bin", by.y = "exp_bin", 
              all = TRUE, suffixes = c("_obs", "_exp"))
  
  setnames(dt, "obs_bin", "bin")
  
  # change any NA from merging to 0
  dt[is.na(N_obs), N_obs:= 0]
  dt[is.na(N_exp), N_exp:= 0]
  
  # add interval boundaries to table
  dt[, `:=`(lb = breaks[-length(breaks)], ub = breaks[-1])]
  
    # THE FOLLOWING WHILE LOOP ISN'T NECESSARY WITH EQUIPROBABLE BINS
  # merge bins (rows) if needed in order to get all(N_exp >= 5) == TRUE

  while(any(dt$N_exp < 5)){
    
    il5 <- max(which(dt$N_exp < 5))
    
    if(il5 == nrow(dt)){
      mrow <- -1
    } else if(il5 == 1) {
      mrow <- +1
    } else {
      if(dt$N_exp[il5 - 1] < dt$N_exp[il5 + 1]){
        mrow <- -1
      } else {
        mrow <- +1
      } 
    }
    dt[il5 + mrow, `:=`(N_obs = dt$N_obs[il5] + dt$N_obs[il5 + mrow],
                        N_exp = dt$N_exp[il5] + dt$N_exp[il5 + mrow])]
    
    # fix bins (just for plotting frequencies)
    if(mrow == -1){
      dt[il5 - 1, bin:= paste0(sub("\\,.*$", "", dt$bin[il5 - 1]),
                               sub("^.*\\,", ",", dt$bin[il5]))]  
    } else {
      dt[il5 + 1, bin:= paste0(sub("\\,.*$", "", dt$bin[il5]),
                               sub("^.*\\,", ",", dt$bin[il5 + 1]))]  
    }
    
    # remove row  
    dt <- dt[-il5, ]
  
  }  # end while loop
  
  # chisquare table values (obs)
  dt[, N_dif:= N_obs - N_exp]
  dt[, chisqt:= (N_dif^2)/N_exp]
  
  # save table for later
  lbins[[i]] <- dt
    
  # compute chi square statistic and save
  chisqv[i] <- dt[, sum(chisqt)]
  
}  

plot(log(crop_area), log(chisqv))
plot(log(crop_niche), log(chisqv))

# chi square values seems to be more associated with niche range than with other factors

lm(log(chisqv) ~ log(crop_niche))$residuals %>% `[`(order(.))



for(i in 1:length(crops)){
  dt <- lbins[[i]] 
  plot(dt$lb, dt$N_dif/dt$N_exp, type = "l", 
       xlab = "bin lower limit", ylab = "(obs - exp)/exp", main = crops[[i]])
  abline(h = 0, lty = 3, col = "red")
}



