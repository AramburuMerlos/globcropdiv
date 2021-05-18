library(magrittr)
library(data.table)
library(terra)
library(Recocrop)

if(system('hostname', TRUE) == "ESP-RH-9891"){
  setwd("D:/globcropdiv/")
} else if(system('hostname', TRUE) == "LAPTOP-ST129J47"){ 
  setwd("G:/My Drive/globcropdiv/")
} # else if { ... 

# Functions ------
# cumulative frequencies 
cf <- function(x) cumsum(x/sum(x))

# Data Prep ------
# Crop Species with Ecocrop Parameters 
dpars <- fread("AuxData/CropSpecies.csv") 

# input data directories
wcdir <- "InData/WorldClim/2.1/wc5min"
sdir <- "InData/SoilGrids"

# Prepare data table
d <- "InData/TotalCropland.tif" %>%
  rast() %>%
  values() %>% 
  is.na() %>% `!` %>% 
  which() %>% 
  data.table(cell = .)

# Crop abundance
abun <- "InData/CropAbundance/*.tif" %>% Sys.glob() %>% rast() 
crops <- names(abun)
d[, (crops):= extract(abun, cell)]

# Crop abundance data quality (for Monfreda crops)
dq <- "InData/CropAbundance/DataQualityIndex/*.tif" %>% Sys.glob() %>% rast() 
dqc <- names(dq)
d[, (paste0(dqc,"_dq")):= extract(dq, cell)]

area_before <- colSums(d[, ..crops])

# For Monfreda crops, keep only high quality data 
# that is, from districts with enviromental CV < 0.5. 
for(j in names(dq)){
 j_dq <- paste0(j, "_dq")
 set(d, i = which(is.na(d[[j_dq]])), j_dq, 1) # change NA data qualities to 1. 
 set(d, j = j, value = ifelse(d[[j_dq]] > 0.5, 0, d[[j]]))
}

area <- colSums(d[, ..crops], na.rm = T)
 
# Static Variables ##########################################################

## pH -----------------------
r_ph <- rast(file.path(sdir, "phh2o/phh2o_0-15cm_mean_5min.tif"))
d[, ph:= extract(r_ph, cell)]
d[, ph:= ph/10]

# cumulative frequency of crop shares 
setorder(d, ph)
pH_cf <- paste0(crops, ".pHcf")
d[, (pH_cf):= lapply(.SD, cf), .SDcols = crops]

# Parameters to be calibrated (only if a more extreme value is detected)
pars <- c("PHMIN", "PHOPMN", "PHOPMX", "PHMAX")
# quantiles used to calibrate each parameter
quan <- c(0.025, 0.25, 0.75, 0.975)
names(quan) <- pars

# calibrate pH
for(i in 1:length(crops)){
  for(j in pars){
    r <- which(dpars[,crop] == crops[i])
    q_pH <- (d[[pH_cf[i]]] - quan[j]) %>%
      abs() %>%
      which.min() %>%
      `[`(d$ph, .)
    if(j %in% c("PHMIN", "PHOPMN")){
      p_pH <- min(dpars[[j]][r]) # the lowest of all spp in that crop group
      if(q_pH < p_pH){
        tc <- dpars[[j]][r] %in% p_pH
        set(dpars, r[tc], j, q_pH)
      }
    } else {
      p_pH <- max(dpars[[j]][r]) # the greatest of all spp in that crop group
      if(q_pH > p_pH){
        tc <- dpars[[j]][r] %in% p_pH
        set(dpars, r[tc], j, q_pH)
      }
    }
  }
}

d[, (pH_cf):= NULL]
d[, ph:= NULL]
setorder(d, cell)

## Annual precipitation ---------------
prec <- "prec/*.tif" %>% file.path(wcdir, .) %>% Sys.glob() %>% rast()
annprec <- app(prec, fun = sum, na.rm = T)

d[, ap:= extract(annprec, cell)]

# add irrigation data 
irrigation <- rast("InData/AQUASTAT/gmia_v5_aei_pct.asc")
d[, irr:= extract(irrigation, cell)]

# set to NA places with more than 10% of irrigation
d[irr > 10, ap:= NA]

# cumulative frequency of crop shares (avoiding NA cells)
setorder(d, ap, na.last = TRUE)
ap_cf <- paste0(crops, ".apcf")
d[!is.na(ap), (ap_cf):= lapply(.SD, cf), .SDcols = crops]


# Parameters to be calibrated (only if a more extreme value is detected)
pars <- c("RMIN", "ROPMN", "ROPMX", "RMAX")
# quantiles used to calibrate each parameter
quan <- c(0.025, 0.25, 0.75, 0.975)
names(quan) <- pars

# calibrate Annual Precipitation
for(i in 1:length(crops)){
  for(j in pars){
    r <- which(dpars[,crop] == crops[i])
    q_ap <- (d[[ap_cf[i]]] - quan[j]) %>%
      abs() %>%
      which.min() %>%
      `[`(d$ap, .)
    if(j %in% c("RMIN", "ROPMN")){
      p_ap <- min(dpars[[j]][r]) # the lowest of all spp in that crop group
      if(q_ap < p_ap){
        tc <- dpars[[j]][r] %in% p_ap
        set(dpars, r[tc], j, q_ap)
      }
    } else {
      p_ap <- max(dpars[[j]][r]) # the greatest of all spp in that crop group
      if(q_ap > p_ap){
        tc <- dpars[[j]][r] %in% p_ap
        set(dpars, r[tc], j, q_ap)
      }
    }
  }
}

d[, (ap_cf):= NULL]
d[, ap:= NULL]
d[, irr:= NULL]
gc(reset = T)
rm(list = ls()[!ls() %in% c("wcdir", "crops", "dpars", "d", "area")])

setorder(d, cell)


# Dynamic Variables ###########################################################
## Parameters --------
# change integer columns to numeric in dpars
intcols <- names(dpars)[sapply(dpars, class) == "integer"]
dpars[, (intcols):= lapply(.SD, as.numeric), .SDcols = intcols]

tpar <- c("TMIN","TOPMN","TOPMX","TMAX")
rpar <- c("RMIN","ROPMN","ROPMX","RMAX")
mrpar <- paste0(rpar, "_M")

# Duration
dpars[, duration:= GMIN + pmin(30, GMAX - GMIN, na.rm=TRUE)]
dpars[, duration:= 15 * round(duration / 15)]

# Monthly rainfall 
# (same method as in Ecocrop, but different starting values)
dpars[, diva:= (GMAX + GMIN) / 2]
dpars[, div1:= ifelse(GMAX > GMIN, diva + 30, diva)]
dpars[, div2:= ifelse(GMAX > GMIN, diva - 30, diva)]
dpars[, RMIN_M:= round(RMIN/(div1/30))]
dpars[, ROPMN_M:= round(ROPMN/(((div1 + diva)/2)/30))]
dpars[, ROPMX_M:= round(ROPMX/(((div2 + diva)/2)/30))]
dpars[, RMAX_M:= round(RMAX/(div2/30))]
dpars[, c("diva", "div1", "div2"):= NULL]


## Upload Variables ------------
# Minimum (Killing) Temperature
r_tmin <- "tmin/*.tif" %>% file.path(wcdir, .) %>% Sys.glob() %>% rast()
tmin <- paste0("tmin_", 1:12) 
d[, (tmin):= extract(r_tmin, cell)]

# Monthly Avg temp
r_tavg <- "tavg/*.tif" %>% file.path(wcdir, .) %>% Sys.glob() %>% rast()
tavg <- paste0("tavg_", 1:12) 
d[, (tavg):= extract(r_tavg, cell)]

# Precipitation
r_prec <- "prec/*.tif" %>% file.path(wcdir, .) %>% Sys.glob() %>% rast()
prec <- paste0("prec_", 1:12) 
d[, (prec):= extract(r_prec, cell)]

# irrigation
r_irri <- rast("InData/AQUASTAT/gmia_v5_aei_pct.asc")
d[, irri:= extract(r_irri, cell)]

## Crop Calibration -------------
# for each crop
# compute scores using model, if prop(scores == 0) > 0.05, change pars 

# parameters that will be calibrated
pcal <- c("KTMP", "TMIN", "TMAX", "RMIN_M", "RMAX_M")

for(i in 1:length(crops)){
  
  ### Default parameters ------
  spp <- which(dpars$crop == crops[i])
  score_spp <- paste0("score_", spp) 
  
  rcci <- d[, .I[d[[crops[i]]] > 0 & irri < 10]] # rainfed cells with crop i
  icci <- d[, .I[d[[crops[i]]] > 0 & irri >= 10]] # irrigated cells with crop i
  
  # run for each species within crop category
  for(j in 1:length(spp)){
    # crop species parameters
    crop_pars <- cbind(
      duration = c(dpars[spp[j], duration], NA, NA, NA),
      ktmp = dpars[spp[j], KTMP] + c(-1, +1, NA, NA),
      tavg = unlist(dpars[spp[j], ..tpar], use.names = FALSE),
      prec = unlist(dpars[spp[j], ..mrpar], use.names = FALSE)
    )

    foo <- list(name = dpars[spp[j], NAME], 
                parameters = crop_pars)
    # rainfed 
    m <- ecocrop(foo)
    control(m, get_max=TRUE)
    
    dynamicPredictors(m) <- cbind(
      ktmp = as.vector(t(d[rcci, ..tmin])),
      tavg = as.vector(t(d[rcci, ..tavg])),
      prec = as.vector(t(d[rcci, ..prec]))
    ) 
    d[rcci, (score_spp[j]):= run(m)]
    
    # irrigated
    m <- ecocrop(foo)
    control(m, get_max=TRUE)
    dynamicPredictors(m) <- cbind(
      ktmp = as.vector(t(d[icci, ..tmin])),
      tavg = as.vector(t(d[icci, ..tavg]))
    ) 
    d[icci, (score_spp[j]):= run(m)]
  }
  
  d[, score:= do.call(pmax, .SD), .SDcols = score_spp]
  d[, (score_spp):= NULL]
  # omission error
  oe <- d[score == 0, sum(.SD), .SDcols = crops[i]] / area[crops[i]]
  
  ## Change parameters  ---------
  while(oe > 0.05){
    d0 <- d[score == 0] 
    par_spp <- paste0("score_", rep(pcal, each = length(spp)), "_sp", spp)
    oa <- rep(NA, length(par_spp))   # vector to store omitted areas
    names(oa) <- gsub("score_", "", par_spp)
    # find the species-parameter combination that produces the greatest change
    for(p in 1:length(pcal)){
      for(j in 1:length(spp)){
        pj <- (p-1) * length(spp) + j
        crop_pars <- cbind(
          duration = c(dpars[spp[j], duration], NA, NA, NA),
          ktmp = dpars[spp[j], KTMP] + c(-1, +1, NA, NA),
          tavg = unlist(dpars[spp[j], ..tpar], use.names = FALSE),
          prec = unlist(dpars[spp[j], ..mrpar], use.names = FALSE)
        )
        # change parameter
        if(pcal[p]=="KTMP") {
          crop_pars[1:2,"ktmp"] <- crop_pars[1:2,"ktmp"] - 1
        }
        if(pcal[p]=="TMIN") {
          crop_pars[1,"tavg"] <- max(crop_pars[1,"tavg"] - 1, 0)
        }
        if(pcal[p]=="TMAX") {
          crop_pars[4,"tavg"] <- crop_pars[4,"tavg"] + 1
        }
        if(pcal[p]=="RMIN_M") {
          crop_pars[1,"prec"] <- max(crop_pars[1,"prec"] - 5, 0)
        }
        if(pcal[p]=="RMAX_M") {
          crop_pars[4,"prec"] <- crop_pars[4,"prec"] + 50
        } 
          
        foo <- list(name = dpars[spp[j], NAME], 
                    parameters = crop_pars)
        # rainfed 
        if(nrow(d0[irri < 10]) > 0){
          m <- ecocrop(foo)
          control(m, get_max=TRUE)
          
          dynamicPredictors(m) <- cbind(
            ktmp = as.vector(t(d0[irri < 10, ..tmin])),
            tavg = as.vector(t(d0[irri < 10, ..tavg])),
            prec = as.vector(t(d0[irri < 10, ..prec]))
          ) 
          d0[irri < 10,(par_spp[pj]):= run(m)]  
        }
        
        # irrigated
        if(nrow(d0[irri > 10]) > 0){
          m <- ecocrop(foo)
          control(m, get_max=TRUE)
          dynamicPredictors(m) <- cbind(
            ktmp = as.vector(t(d0[irri >= 10, ..tmin])),
            tavg = as.vector(t(d0[irri >= 10, ..tavg]))
          ) 
          d0[irri >= 10, (par_spp[pj]):= run(m)]
        }
        # omitted area for species parameter combination
        oa[pj] <- sum(d0[[crops[i]]][d0[[par_spp[pj]]] == 0])
      }
    }
    # which species-param combinations resulted in less omitted area
    wpj <- which.min(oa) %>% names()

    # update score in main data.table and omission error
    d[score == 0, score:= d0[[paste0("score_", wpj)]]]
    oe_u <- oa[wpj]/area[crops[i]]
    
    # if improvement, update omission error and parameter, else, break loop
    if(oe_u < oe){
      oe <- oe_u
      wsp <- sub(".*sp", "", wpj) %>% as.numeric() # winner species 
      wp <- sub("_sp.*$", "", wpj)                 # winner parameter  
      if(wp == "KTMP") dpars[wsp, KTMP:= KTMP - 1]
      if(wp == "TMIN") dpars[wsp, TMIN:= max(TMIN - 1, 0)]
      if(wp == "TMAX") dpars[wsp, TMAX:= TMAX + 1]
      if(wp == "RMIN_M") dpars[wsp, RMIN_M:= max(RMIN_M - 5, 0)]
      if(wp == "RMAX_M") dpars[wsp, RMAX_M:= RMAX_M + 50]
    } else {
      warning(paste("no solution for", crops[i]))
      break
    }
  }
  print(paste0(i,"/", length(crops)))
}

fwrite(dpars, "AuxData/CalibEcoPars.csv") 
