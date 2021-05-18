library(magrittr)
library(data.table)
library(terra)
library(ranger)

if(system('hostname', TRUE) == "ESP-RH-9891"){
  setwd("D:/globcropdiv/")
} else if(system('hostname', TRUE) == "LAPTOP-ST129J47"){ 
  setwd("G:/My Drive/globcropdiv/")
} # else if { ... 


# FUNCTIONS ###########################################################

frmse <- function(obs, pred) sqrt(mean((obs-pred)^2, na.rm = T))


# DATA PREP ###########################################################

## Diversity data (y) ----------
Da_r <- rast("OutData/Da.tif")
Dp_int_r <- rast("OutData/Dp_int.tif")
Dp_sdm_r <- rast("OutData/Dp_sdm.tif")
Dp_eco_r <- rast("OutData/Dp_eco.tif")

d <- rast("InData/TotalCropland.tif") %>% 
  values() %>% {data.table(cell = which(!is.na(.)))}

d[, Da:= extract(Da_r, cell)]
d[, Dp_int:= extract(Dp_int_r, cell)]
d[, Dp_sdm:= extract(Dp_sdm_r, cell)]
d[, Dp_eco:= extract(Dp_eco_r, cell)]

y <- c("Da", "Dp_int", "Dp_sdm", "Dp_eco")

## Predictors --------------
preds_r <- fread("AuxData/SelectedPreds.csv")$file_name %>% 
  file.path("InData", .) %>% 
  Sys.glob() %>% 
  rast()
names(preds_r) <- preds <- fread("AuxData/SelectedPreds.csv")$short_name
d[, (preds):= extract(preds_r, cell)]

r_field <- rast(Sys.glob("InData/field_size*/field_size.tif"))
d[, field_size:= extract(r_field, cell)]
preds <- c(preds, "field_size")


# TUNE MODELS #################################################################
# 5% of the data will be used for model tuning, sampled regularly 
tu_ss = 0.05
itu <- which(d$cell%%(1/tu_ss) == 0)
dtu <- d[itu,]

# tuning parameters
tp <- expand.grid(
  mtry = seq(3, 15, 3),
  min.node.size = c(1, seq(20, 100, 20)),
  Da_oob_rmse = NA_real_, 
  Dp_int_oob_rmse = NA_real_,
  Dp_sdm_oob_rmse = NA_real_, 
  Dp_eco_oob_rmse = NA_real_
)
setDT(tp)


set.seed(0)
for(i in 1:nrow(tp)){
  for(j in y){
    rf <- ranger(x = dtu[!is.na(dtu[[j]]), ..preds], 
                 y = dtu[!is.na(dtu[[j]])][[j]],
                 num.trees = 1000,
                 mtry = tp$mtry[i],
                 min.node.size = tp$min.node.size[i])
    tp[i, paste0(j, "_oob_rmse"):= rf$prediction.error]
    print(paste0("RF: ", i, "/", nrow(tp), " - ", j))
  }
}

fwrite(tp, "OutData/RF_tuning.csv")


