library(magrittr)
library(data.table)
library(terra)

if(!getwd() %like% "globcropdiv$") warning("See 0000_wd.R")

source("Functions/allocate.R")

# Data prep -------

d <- "InData/TotalCropland.tiff" %>%
  rast() %>%
  values() %>% 
  is.na() %>% `!` %>% 
  which() %>% 
  data.table(cell = .)

# Suitability
suit <- "OutData/Int_Suit/*.tif" %>% Sys.glob() %>% rast()
crops <- names(suit)
cs <- paste0(crops, ".suit")
names(suit) <- cs

# Abundance
abun <- fread("AuxData/CropAbundanceSource.csv") %>%
  merge(data.table(crop = crops)) %$% file %>% 
  file.path("D:", .) %>% rast()

ca <- paste0(crops, ".abun")
names(abun) <- ca

# read data
d[, (ca):= extract(abun, cell)]
d[, (cs):= extract(suit, cell)]
 
# change any NA to 0
for(j in names(d)) set(d, which(is.na(d[[j]])), j, 0)

# remove cells unsuitable for all crops (there shouldn't be any, but there are)
d <- d[rowSums(d[,..cs] > 0) > 0,]

## parameters -----
#  starting point for allocation of cell cropland to crop i for initial rounds
start.alloc <- 0.05

# max allocation(frac):suitability(frac) ratio: 
## the proportion of allocated land to crop j can't be more than as_rat times
## its normalized suitability.
as_rat <- 100

d <- allocate(d, ca, cs, start.alloc, as_rat)

fwrite(d, "OutData/allocated_v01.csv")

# Cross Entropy function
cef <- function(d, suitcols, sharecols){
  nc <- length(suitcols)
  c1 <- c2 <- rep(0, nc)
  for(j in 1:nc){
    suit <- d[[suitcols[j]]]
    share <- d[[sharecols[j]]]
    rm0 <- suit == 0 | share == 0
    c1[j] <- sum(share[!rm0] * log(share[!rm0]))
    c2[j] <- sum(share[!rm0] * log(suit[!rm0]))
  }
  return(sum(c1) - sum(c2))
}
names(d)
summary(d$ra)
sum(d$maize.allo)
sum(d$maize.abun)
cor(d$maizefor.allo / d$cl, d$maize.abun)

aa <- paste0(crops, ".allo")
plot(colSums(d[,..abun_columns]), colSums(d[,..allo_columns]))
abline(a=0, b=1)


# crop fraction
cf <- gsub(".suit", ".frac", cs)
d[, (cf):= lapply(.SD, function(x) x/cl), .SDcols = ]

cef(d, cs, cf)

d

# add reallocation algorithm trying to minimize cross entropy


# OLD R#########################################################################


# d[1, maize.suit:=0]

# set area to 0 when suit is 0 (suit is never 0 in Maxent?)
for(j in 1:ncr) set(d, which(d[[crops.suit[j]]] == 0), crops.area[j], 0)

# compute cell shares of crop area
crops.share <- gsub(".abun", ".share", crops.abun)

for(j in 1:ncr) d[, (crops.share[j]):= .SD/a[j], .SDcols = crops.area[j]]


# column to keep remaining cropland per cell
d[, rc:= cl - Reduce(`+`, .SD), .SDcols = crops.area]


# Cross Entropy function
cef <- function(d, suitcols, sharecols){
  nc <- length(suitcols)
  c1 <- c2 <- rep(0, nc)
  for(j in 1:nc){
    suit <- d[[suitcols[j]]]
    share <- d[[sharecols[j]]]
    rm0 <- suit == 0 | share == 0
    c1[j] <- sum(share[!rm0] * log(share[!rm0]))
    c2[j] <- sum(share[!rm0] * log(suit[!rm0]))
  }
  return(sum(c1) - sum(c2))
}


ce.bf <- cef(d, crops.suit, crops.share)
ce.af = ce.bf


maxiter <- 100

prev.d <- data.table::copy(d)
prev.d
iter = 0

while(iter < maxiter){
  if(ce.af < ce.bf){
    iter = 1
    ce.bf <- ce.af
  } else {
    iter = iter + 1
    d <- data.table::copy(prev.d)
  }
  # relocation algorithm -----
  # remove 10/iter percent of the crop area for all crops
  d[, (crops.area):= .SD * (1 - 0.1/iter), .SDcols = (crops.area)]
  d[, rc:= cl - Reduce(`+`, .SD), .SDcols = crops.area]
  
  # crop area to be reallocated
  rca <- a - d[, lapply(.SD, sum), .SDcol = (crops.area)]
  
  # relocate crop area to greatest suitability, in random order
  rand.ord <- sample(1:ncr)
  
  for(j in rand.ord){
    setorderv(d, cols = crops.suit[j], order = -1L)
    # area to fill based on suitability
    d[, cum.area := cumsum(rc * .SD), .SDcols = crops.suit[j]]
    # enough cropland to fill cells until istop - 1
    istop <- Position(function(x) x >= rca[[j]], d$cum.area)
    if(is.na(istop)) istop <- nrow(d) 
    d[0:(istop - 1), (crops.area[j]):= .SD + rc, .SDcols = crops.area[j]]
    # remainder
    rcaj <- max(c(a[j] - d[, sum(.SD), .SDcol = crops.area[j]],0))
    d[istop, (crops.area[j]):= .SD + rcaj, .SDcols = crops.area[j]]
    
    #update remaining cropland per cell
    d[, rc:= cl - Reduce(`+`, .SD), .SDcols = crops.area]
  }
  
  # update shares
  for(j in 1:ncr) d[, (crops.share[j]):= .SD/a[j], .SDcols = crops.area[j]]
  ce.af <- cef(d, crops.suit, crops.share)
}



d











# number of ha to be assigned per loop (can be progressive)
xs <- 2^(5:-4)

# x ha at a time assigned to  a single cell takes forever. Maybe assign more cells simultaneously (z % of available cell  and y % of remaining crop area?q)



for(x in xs){
  while(any(d$rc > x)){
    # sample a crop with more than x ha to be assigned
    j <- 
      d[, lapply(.SD, sum), .SDcols = crops.area] %>%
      unlist(use.names = F) %>%
      `-`(ca, .) %>%{
        sample.vec(which(. > x), size = 1L) #, prob = .)   
      }  
    # sample a cell with free cropland using suitability as prob
    i <- 
      d[, .I[rc > x]] %>%  
      sample.vec(size = 1L, prob = d[[crops.suit[j]]][.])
    
    # set x cropland area of cell i to crop j
    set(d, i, crops.area[j], d[[crops.area[j]]][i] + x)
    
    # compute remaining cropland
    set(d, i, "rc", d$rc[i] - x)
  }
}


d[, sum(cl)]/32 / d[, sum(cl)-sum(rc)]/32


d[, ts:= Reduce(`+`, .SD), .SDcols = crops.suit]




# suit x tc per cell
crops.sxtc <- paste0(crops, ".sxtc")
d[, (crops.sxtc):= lapply(.SD, function(x) x * tc), .SDcols =  crops.suit]







# total cropland
totcl <- rast("InData/TotalCropland.tiff")
d <- data.table(cell = 1:ncell(totcl), cl = values(totcl)[,1])
d <- d[!is.na(cl),]


d[, (ns):= extract(suit, cell)] # might be faster if saved in SSD
d <- d[complete.cases(d), ]

d[, (na):= extract(abun, cell)]

d[, lapply(.SD, sum), .SDcols = (nra)]

setDT(rarea, keep.rownames = "crops")
setnames(rarea, "sum", "real_area")

area <- scarea[rarea, on = "crops"]
area[, ratio:= sim_area/real_area]


# sum of all suitability values per cell
d[, ts:= Reduce('+', .SD), .SDcols = (ns)]

# simulated abundance per cell
nsa <- gsub("suit", "sabun", ns)
d[, (nsa):= lapply(.SD, function(x) x * cl), .SDcols = (ns)]


# get total (simulated) crop area 
scarea <- global(sabun, "sum", na.rm = T)
setDT(scarea, keep.rownames = "crops")
setnames(scarea, "sum", "sim_area")

# over predicted crops
op_crops <- area[ratio > 1, crops]
# under predicted crops
un_crops <- area[ratio < 1, crops]

# reduce area for overpredicted crops proportionally across cells
fabun <- sabun 

for(crop in op_crops){
  fabun[[crop]] <- sabun[[crop]] / area[crops == crop, ratio]
}  





diver <- function(x){
  p <- x / sum(x)
  H <- -sum(p * log(p))
  return(exp(H))
}

area[, lapply(.SD,diver), .SDcols = c("sim_area", "real_area")]


# the last step should "free" cropland 

# for crops which estimated crop area doesn't reach their total crop area
# increase their abundance in cells with "free" cropland
# how? start with crops whith greatest demand for cropland?
# assign free cropland to that crop based on suitability to fulfill total crop area
# recalculate free cropland, move to next crop


