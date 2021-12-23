library(magrittr)
library(terra)
library(data.table)

if(system('hostname', TRUE) == "ESP-RH-9891"){
  setwd("D:/globcropdiv/")
} else if(system('hostname', TRUE) == "LAPTOP-ST129J47"){ 
  setwd("G:/My Drive/globcropdiv/")
} # else if { ... 

countries <- rast("OutData/projected/rcountries_proj.tif")

# diversity function 
fd <- function(x) exp(-sum(x * log(x), na.rm = T))

cl <- rast("OutData/projected/CroplandProp.tif")

cells <- which(!is.na(values(cl)))

# Local Average country diversity (Dalpha) ###########
divs <- c("act_D", "att_D_eco", "att_D_sdm")

rD <- rast(paste0("OutData/", divs, ".tif"))



# extract values
d <- data.table(cell = cells, 
                pcl = extract(cl, cells)[,1], 
                country = extract(countries, cells)[,1])
d[, (divs):= extract(rD, cell)]
d[, tcl:= pcl * prod(res(cl))/1e4]

# Average attainable diversity
d[, att_D_avg:= (att_D_eco + att_D_sdm)/2]

# add contintent
codes <- geodata::country_codes()
setDT(codes)
codes <- codes[, .SD, .SDcols = c("ISO3", "continent")]
setnames(codes, "ISO3", "country")
codes[continent %like% "America", continent:= "America"]

d <- codes[d, on = c("country")]


# compute D alpha
fda <- function(D, tcl)  exp(sum(log(D) * tcl/sum(tcl, na.rm = T), na.rm = T))

da <- d[, lapply(.SD, fda, tcl = tcl), .SDcols = c('act_D', 'att_D_avg'), by = continent]
da <- da[!is.na(continent),]

da[, Dg:= (att_D_avg - act_D)/att_D_avg * 100]
da
