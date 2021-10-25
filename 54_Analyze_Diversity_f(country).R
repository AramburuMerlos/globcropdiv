library(terra)
library(data.table)
library(magrittr)

if(system('hostname', TRUE) == "ESP-RH-9891"){
  setwd("D:/globcropdiv/")
} else if(system('hostname', TRUE) == "LAPTOP-ST129J47"){ 
  setwd("G:/My Drive/globcropdiv/")
} # else if { ... 

countries <- rast("InData/countries/rcountries.tif")

# diversity function 
fd <- function(x) exp(-sum(x * log(x), na.rm = T))

totcl <- rast("InData/TotalCropland.tif")

# Total diversity per countries (Dgamma) ##########################

## Actual diversities ########

actual <- "InData/CropAbundance/*.tif" %>%  Sys.glob() %>%  rast()
crops <- names(actual)

# extract values, long 
cells <- values(totcl) %>% {which(!is.na(.))}

l <- vector(mode = "list", length = length(crops))

for(i in 1:length(crops)){
  dt <- data.table(crop = crops[i], 
                   area = extract(actual[[i]], cells)[,1], 
                   country = extract(countries, cells)[,1])
  dt <- dt[!is.na(area),]
  dt <- dt[area > 0,]
  l[[i]] <- dt
  rm(dt)
}

d <- rbindlist(l)
rm(l); gc(reset = T)

d <- d[, .(area = sum(area)), by = .(country, crop)]

d[, prop:= area/sum(area), by = .(country)]

dd <- d[, .(tot.area = sum(area), act_Dg = fd(prop)), by = .(country)]
dd <- dd[!is.na(country), ]

rm(d); gc(reset = T)

## Attainable diversity ############

### Ecocrop ##############
alloc_eco <- fread("OutData/allocated_eco.csv")
d <- melt(alloc_eco, id.vars = "cell", measure.vars = patterns("^a."),
          variable.name = "crop", value.name = "area")
d <- d[area > 0, ]

d[, country:= extract(countries, cell)[,1]]
d <- d[, .(area = sum(area)), by = .(country, crop)]
d[, prop:= area/sum(area), by = .(country)]

d <- d[, .(att_Dg_eco = fd(prop)), by = .(country)]
dd <- d[dd, on = "country"]

rm(d, alloc_eco); gc(reset = T)


### SDM ##############
alloc_sdm <- fread("OutData/allocated_sdm.csv")
d <- melt(alloc_sdm, id.vars = "cell", measure.vars = patterns("^a."),
          variable.name = "crop", value.name = "area")
d <- d[area > 0, ]

d[, country:= extract(countries, cell)[,1]]
d <- d[, .(area = sum(area)), by = .(country, crop)]
d[, prop:= area/sum(area), by = .(country)]

d <- d[, .(att_Dg_sdm = fd(prop)), by = .(country)]
dd <- d[dd, on = "country"]

rm(d, alloc_sdm); gc(reset = T)


## Potential Diversity ########################

### Ecocrop ############################
eco_suit <- "OutData/Ecocrop/*.tif" %>%  Sys.glob() %>%  rast()
all.equal(names(eco_suit), crops)

# extract values, long 
l <- vector(mode = "list", length = length(crops))

for(i in 1:length(crops)){
  dt <- data.table(crop = crops[i],
                   cell = cells,
                   tcl = extract(totcl, cells)[,1],
                   suit = extract(eco_suit[[i]], cells)[,1], 
                   country = extract(countries, cells)[,1])
  dt <- dt[!is.na(suit),]
  dt <- dt[suit > 0,]
  l[[i]] <- dt
  rm(dt)
}

d <- rbindlist(l)
rm(l); gc(reset = T)

# normalize suitabilities
d[, suit:= suit/sum(suit, na.rm = T), by = .(cell)]

# compute areas
d[, area:= suit * tcl]

# total crop areas by country 
d <- d[, .(area = sum(area)), by = .(country, crop)]

d[, prop:= area/sum(area), by = .(country)]

d <- d[, .(pot_Dg_eco = fd(prop)), by = .(country)]
d <- d[!is.na(country),]

dd <- d[dd, on = "country"]


### SDM ############################
sdm_suit <- "OutData/SDM/*_AVG.tif" %>%  Sys.glob() %>%  rast()
all.equal(names(sdm_suit), crops)

# extract values, long 
l <- vector(mode = "list", length = length(crops))

for(i in 1:length(crops)){
  dt <- data.table(crop = crops[i],
                   cell = cells,
                   tcl = extract(totcl, cells)[,1],
                   suit = extract(sdm_suit[[i]], cells)[,1], 
                   country = extract(countries, cells)[,1])
  dt <- dt[!is.na(suit),]
  dt <- dt[suit > 0,]
  l[[i]] <- dt
  rm(dt)
}

d <- rbindlist(l)
rm(l); gc(reset = T)

# normalize suitabilities
d[, suit:= suit/sum(suit, na.rm = T), by = .(cell)]

# compute areas
d[, area:= suit * tcl]

# total crop areas by country 
d <- d[, .(area = sum(area)), by = .(country, crop)]

d[, prop:= area/sum(area), by = .(country)]

d <- d[, .(pot_Dg_sdm = fd(prop)), by = .(country)]
d <- d[!is.na(country),]

dd <- d[dd, on = "country"]



rm(d, eco_suit, sdm_suit, actual); gc(reset = T)





# Local Average country diversity (Dalpha) ###########
divs <- c("act_D", "att_D_eco", "att_D_sdm", "pot_D_eco", "pot_D_sdm", "pot_D_eco_sp")

rD <- rast(paste0("OutData/", divs, ".tif"))

# extract values
d <- data.table(cell = cells, 
                tcl = extract(totcl, cells)[,1], 
                country = extract(countries, cells)[,1])
d[, (divs):= extract(rD, cell)]

# compute D alpha
fda <- function(D, tcl)  exp(sum(log(D) * tcl/sum(tcl, na.rm = T), na.rm = T))

da <- d[, lapply(.SD, fda, tcl = tcl), .SDcols = divs, by = country]
da <- da[!is.na(country),]

adivs <- gsub("D", "Da", divs)

setnames(da, divs, adivs)

dd <- da[dd, on = "country", all = T]

rm(d, da); gc(reset = T)


fwrite(dd, "OutData/D_by_country.csv")
