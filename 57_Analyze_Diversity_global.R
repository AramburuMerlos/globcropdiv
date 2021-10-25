library(data.table)
library(terra)

if(system('hostname', TRUE) == "ESP-RH-9891"){
  setwd("D:/globcropdiv/")
} else if(system('hostname', TRUE) == "LAPTOP-ST129J47"){ 
  setwd("G:/My Drive/globcropdiv/")
} # else if { ... 

totcl <- values(rast("InData/TotalCropland.tif"))
cells <- which(!is.na(totcl))
totcl <- totcl[cells]
act_D <- extract(rast("OutData/act_D.tif"), cells)[,1]

cD <- exp(sum(totcl/(sum(totcl)) * log(act_D)))


att_D <- (extract(rast("OutData/att_D_eco.tif"), cells)[,1] + 
            extract(rast("OutData/att_D_eco.tif"), cells)[,1])/2

att_D[is.na(att_D)] <- 1

aD <- exp(sum(totcl/(sum(totcl)) * log(att_D)))
