# ideas:
# use a model to describe diversity patterns. 
# That is, how much of the variation in diversity (Da, Dp) can be explained by biophysical variables.
# Other factors to considered: field size, crop composition (how dominant the crops within that cell are), crop suitability
# analyze variable importance with random forest or something like that
library(terra)
library(data.table)
library(magrittr)
library(xgboost)

if(!grepl("globcropdiv$", getwd())){
  if (system('hostname', TRUE) %in% c("ESP-RH-9891", "LAPTOP-ST129J47")) { 
    setwd("G:/My Drive/globcropdiv/")
  } # else if { ... 
}


Da_r <- rast("OutData/Da.tif")
Dp_int_r <- rast("OutData/Dp_int.tif")

d <- rast("InData/TotalCropland.tif") %>% 
  values() %>% 
  {data.table(cell = which(!is.na(.)), tcl = .[!is.na(.)])}

d[, Da:= extract(Da_r, cell)]
d[, Dp_int:= extract(Dp_int_r, cell)]


d[, lat:= yFromCell(Da_r, cell)]
lat_avg <- d[, .(Da = sum(Da * tcl/sum(tcl)), 
                 Dp_int = sum(Dp_int * tcl/sum(tcl), na.rm = T)), 
             by = lat]

plot(lat_avg[, c(3,1)])
plot(lat_avg[, 2:1])


