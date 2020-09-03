library(here)
library(terra)

# for crop species with multiples entries in the EcoCrop Data base, 
# keep the highest suitability value for the species. 
crops <- read.csv(here("AuxData/selsp.csv"), stringsAsFactors = F) 
dup_sid <- unique(crops$SID[duplicated(crops$SID)])

fn <- Sys.glob("D:/Dp_world/LowRes_10min/CropSuit/*.tif")
suit <- rast(fn)

for(i in 1:length(dup_sid)){
  ri <- which(crops$SID == dup_sid[i]) # crops rows with species i
  # max values for the species 
  smax <- app(suit[[ri]], max)
  # remove all old species files
  file.remove(fn[ri])
  # save species max suit 
  sid <- paste0("ID", formatC(dup_sid[i], width = 4, flag = "0"))
  fn_i <- file.path("D:/Dp_world/LowRes_10min/CropSuit", paste0(sid, ".tif")) 
  writeRaster(smax, filename = fn_i, 
              overwrite = T, 
              wopt = list(names = sid, filetype = "GTiff"))
}


