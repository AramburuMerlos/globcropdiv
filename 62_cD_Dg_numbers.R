library(data.table)
library(terra)


if(system('hostname', TRUE) %in% c("ESP-RH-9891", "LAPTOP-ST129J47")){
  setwd("D:/globcropdiv/")
} # else if { ... 




pcl <- rast("OutData/projected/CroplandProp.tif")
totcl <- pcl * prod(res(pcl))/1e4

act_D <- rast("OutData/act_D.tif")
att_D_sdm <- rast("OutData/att_D_sdm.tif")
att_D_eco <- rast("OutData/att_D_eco.tif")

att_D_avg <- (att_D_eco + att_D_sdm)/2
Dgap_avg <- (att_D_avg - act_D)/(att_D_avg) * 100

d <- data.table(
  cl = values(totcl)[,1], 
  act = values(act_D)[,1],
  gap = values(Dgap_avg)[,1]
)

d <- d[!is.na(cl)]

d[gap < 0, gap:= 0] 
d[gap > 50, sum(cl)]/sum(d$cl)

d[act < 5, sum(cl)]/sum(d$cl)
d[act < 5 & gap > 75, sum(cl)]/d[act < 5, sum(cl)]
d[act < 5 & gap < 50, sum(cl)]/d[act < 5, sum(cl)]

d[act > 12, sum(cl)]/d[, sum(cl)]
d[act > 12 & gap < 60, sum(cl)]/d[act > 12, sum(cl)]
d[act > 12 & gap < 50, sum(cl)]/d[act > 12, sum(cl)]
