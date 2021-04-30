r <- rast("InData/TotalCropland.tif")
d[, lat:= xFromCell(r, cell)]
div_lat <- d[, .(Dp = mean(Dp), Da = mean(Da)), by = lat]
setorder(div_lat, lat)
plot(div_lat$Da, div_lat$lat, type = 'l', xlim = c(0,50))
lines(div_lat$Dp, div_lat$lat, col = "blue")


# fiel data 
# downloaded from geo wiki app
r_field <- rast(Sys.glob("InData/field_size*/field_size*.img"))

r_field <- expand(r_field, r_tcl)
r_field <- aggregate(r_field, fact = 10, fun = "median", na.rm = TRUE)
compareGeom(r_tcl, r_field)

d[, field_size:= extract(r_field, cell)]
d[field_size == 0, field_size := NA]
d[, Dg:= pmax(1 - Da/Dp, 0)]

plot(r_field)
