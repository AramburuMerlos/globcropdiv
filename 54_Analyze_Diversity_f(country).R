library(magrittr)
library(data.table)
library(terra)

if(system('hostname', TRUE) == "ESP-RH-9891"){
  setwd("D:/globcropdiv/")
} else if(system('hostname', TRUE) == "LAPTOP-ST129J47"){ 
  setwd("G:/My Drive/globcropdiv/")
} # else if { ... 

# diversity function 
fd <- function(x) exp(-sum(x * log(x), na.rm = T))


cl <- rast("OutData/projected/CroplandProp.tif")

# countries -----------
world <- geodata::world(resolution = 1, path = "InData/countries")
world <- project(world, cl)
rasterize(world, cl, "GID_0", 
          filename = "OutData/projected/rcountries_proj.tif",
          overwrite = T, 
          wopt = list(names = "countries", filetype = "GTiff",
                      gdal=c("COMPRESS=Deflate","PREDICTOR=1","ZLEVEL=6")))

rworld <- rast("OutData/projected/rcountries_proj.tif")

# Total diversity per countries (Dgamma) ##########################

# uplaod proportions 
actual <- rast("OutData/projected/ActualCropProp.tif")
pot_eco <- rast("OutData/projected/PotEcoCropProp.tif")
pot_sdm <- rast("OutData/projected/PotSDMCropProp.tif")
att_eco <- rast("OutData/projected/AttEcoCropProp.tif")
att_sdm <- rast("OutData/projected/AttSDMCropProp.tif")

crops <- names(actual)

# create data.table with cell, country and total cropland per cell 
d <- data.table(cell = which(!is.na(values(cl))))
d[, country:= terra::extract(rworld, cell)]
d[, cl_prop:= extract(cl, cell)]
d[, totcl:= cl_prop * (prod(res(cl)))/1e4] # in ha


## Actual diversity ------

# extract values, wide
d[, (crops):= extract(actual, cell)]
d[, (crops):= lapply(.SD, function(x) x*totcl), .SDcols = crops]

# reshape to long
dl <- melt(d, measure.vars = crops, variable.name = "crop", value.name = "area")
dl <- dl[!is.na(area),]
dl <- dl[area > 0,]
dl <- dl[!is.na(country),]

dl <- dl[, .(area = sum(area)), by = .(country, crop)]
dl[, prop:= area/sum(area), by = .(country)]

# compute diversity
dd <- dl[, .(tot.area = sum(area), act_Dg = fd(prop)), by = .(country)]


## Attainable diversity ---------

### Ecocrop -------
# extract values, wide
d[, (crops):= extract(att_eco, cell)]
d[, (crops):= lapply(.SD, function(x) x*totcl), .SDcols = crops]

# reshape to long
dl <- melt(d, measure.vars = crops, variable.name = "crop", value.name = "area")
dl <- dl[!is.na(area),]
dl <- dl[area > 0,]
dl <- dl[!is.na(country),]

dl <- dl[, .(area = sum(area)), by = .(country, crop)]
dl[, prop:= area/sum(area), by = .(country)]

# compute diversity
dl <- dl[, .(tot.area = sum(area), att_Dg_eco = fd(prop)), by = .(country)]

dd <- dl[dd, on = "country"]

gc(reset = T)


### SDM -------
# extract values, wide
d[, (crops):= extract(att_sdm, cell)]
d[, (crops):= lapply(.SD, function(x) x*totcl), .SDcols = crops]

# reshape to long
dl <- melt(d, measure.vars = crops, variable.name = "crop", value.name = "area")
dl <- dl[!is.na(area),]
dl <- dl[area > 0,]
dl <- dl[!is.na(country),]

dl <- dl[, .(area = sum(area)), by = .(country, crop)]
dl[, prop:= area/sum(area), by = .(country)]

# compute diversity
dl <- dl[, .(tot.area = sum(area), att_Dg_sdm = fd(prop)), by = .(country)]

dd <- dl[dd, on = "country"]

gc(reset = T)

## Potential diversity ---------

### Ecocrop -------
# extract values, wide
d[, (crops):= extract(pot_eco, cell)]
d[, (crops):= lapply(.SD, function(x) x*totcl), .SDcols = crops]

# reshape to long
dl <- melt(d, measure.vars = crops, variable.name = "crop", value.name = "area")
dl <- dl[!is.na(area),]
dl <- dl[area > 0,]
dl <- dl[!is.na(country),]

dl <- dl[, .(area = sum(area)), by = .(country, crop)]
dl[, prop:= area/sum(area), by = .(country)]

# compute diversity
dl <- dl[, .(tot.area = sum(area), pot_Dg_eco = fd(prop)), by = .(country)]

dd <- dl[dd, on = "country"]

gc(reset = T)


### SDM -------
# extract values, wide
d[, (crops):= extract(pot_sdm, cell)]
d[, (crops):= lapply(.SD, function(x) x*totcl), .SDcols = crops]

# reshape to long
dl <- melt(d, measure.vars = crops, variable.name = "crop", value.name = "area")
dl <- dl[!is.na(area),]
dl <- dl[area > 0,]
dl <- dl[!is.na(country),]

dl <- dl[, .(area = sum(area)), by = .(country, crop)]
dl[, prop:= area/sum(area), by = .(country)]

# compute diversity
dl <- dl[, .(tot.area = sum(area), pot_Dg_sdm = fd(prop)), by = .(country)]

dd <- dl[dd, on = "country"]

gc(reset = T)




# Local Average country diversity (Dalpha) ###########
divs <- c("act_D", "att_D_eco", "att_D_sdm", "pot_D_eco", "pot_D_sdm")

rD <- rast(paste0("OutData/", divs, ".tif"))

# extract values
d[, (crops):= NULL]
d[, (divs):= extract(rD, cell)]

# compute D alpha
fda <- function(D, tcl)  exp(sum(log(D) * tcl/sum(tcl, na.rm = T), na.rm = T))

da <- d[, lapply(.SD, fda, tcl = totcl), .SDcols = divs, by = country]
da <- da[!is.na(country),]

adivs <- gsub("D", "Da", divs)

setnames(da, divs, adivs)

dd <- da[dd, on = "country", all = T]

rm(d, da); gc(reset = T)

# remove extra total areas
ctrm <- names(dd)[names(dd) %like% "i.tot.area"]
dd[, (ctrm):= NULL]


fwrite(dd, "OutData/D_by_country.csv")
