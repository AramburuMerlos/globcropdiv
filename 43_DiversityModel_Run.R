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


# RUN MODELS ---------------
tp <- fread("OutData/RF_tuning.csv")

## Da ---------------
param_Da <- tp[which.min(Da_oob_rmse),]
Da_cols <- c("Da", preds)

rf_Da <- ranger(data = d[, ..Da_cols], 
                dependent.variable.name = "Da",
                num.trees = 1000,
                mtry = param_Da$mtry,
                min.node.size = param_Da$min.node.size,
                num.threads = 7,
                save.memory = TRUE, 
                importance = "impurity",
                seed = 1234)

saveRDS(rf_Da, "OutData/Da_rf_analysis.RDS")


## Dp_int ---------------
param_Dp_int <- tp[which.min(Dp_int_oob_rmse),]
Dp_int_cols <- c("Dp_int", preds)

rf_Dp_int <- ranger(data = d[!is.na(Dp_int), ..Dp_int_cols], 
                    dependent.variable.name = "Dp_int",
                    num.trees = 1000,
                    mtry = param_Dp_int$mtry,
                    min.node.size = param_Dp_int$min.node.size,
                    num.threads = 7,
                    save.memory = TRUE, 
                    importance = "impurity",
                    seed = 1234)

saveRDS(rf_Dp_int, "OutData/Dp_int_rf_analysis.RDS")


## Dp_sdm ---------------
param_Dp_sdm <- tp[which.min(Dp_sdm_oob_rmse),]
Dp_sdm_cols <- c("Dp_sdm", preds)

rf_Dp_sdm <- ranger(data = d[, ..Dp_sdm_cols], 
                    dependent.variable.name = "Dp_sdm",
                    num.trees = 1000,
                    mtry = param_Dp_sdm$mtry,
                    min.node.size = param_Dp_sdm$min.node.size,
                    num.threads = 7,
                    save.memory = TRUE, 
                    importance = "impurity",
                    seed = 1234)

saveRDS(rf_Dp_sdm, "OutData/Dp_sdm_rf_analysis.RDS")


## Dp_eco ---------------
param_Dp_eco <- tp[which.min(Dp_eco_oob_rmse),]
Dp_eco_cols <- c("Dp_eco", preds)

rf_Dp_eco <- ranger(data = d[!is.na(Dp_eco), ..Dp_eco_cols], 
                    dependent.variable.name = "Dp_eco",
                    num.trees = 1000,
                    mtry = param_Dp_eco$mtry,
                    min.node.size = param_Dp_eco$min.node.size,
                    num.threads = 7,
                    save.memory = TRUE, 
                    importance = "impurity",
                    seed = 1234)

saveRDS(rf_Dp_eco, "OutData/Dp_eco_rf_analysis.RDS")



