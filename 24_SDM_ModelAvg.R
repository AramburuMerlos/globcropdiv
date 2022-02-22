library(magrittr)
library(terra)

if(system('hostname', TRUE) %in% c("ESP-RH-9891", "LAPTOP-ST129J47")){
  setwd("D:/globcropdiv/")
} # else if { ... 


me <- "OutData/SDM/*_ME.tif" %>% Sys.glob() %>% rast()
rf <- "OutData/SDM/*_RF.tif" %>% Sys.glob() %>% rast()
br <- "OutData/SDM/*_BRT.tif" %>% Sys.glob() %>% rast()

all.equal(names(me), names(rf))
all.equal(names(me), names(br))

crops <- names(me)

w <- readRDS("OutData/SDM/CV/Tuning/Model_averaging_weights.RDS")

favg <- function(x,y,z) x * w[1] + y * w[2] + z * w[3]

for(i in 1:length(crops)){
  fn <- file.path("OutData/SDM", paste0(crops[i], "_AVG.tif"))
  wopts <- list(names = crops[i], filetype = "GTiff",
                gdal=c("COMPRESS=Deflate","PREDICTOR=1", "ZLEVEL=6"))
  lapp(c(me[[i]], rf[[i]], br[[i]]), 
       favg, 
       filename = fn, overwrite = T, wopt = wopts)
}
