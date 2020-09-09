library(terra)

# for crop species with multiples entries in the EcoCrop Data base, 
# keep the highest suitability value for the species. 
crops <- read.csv("AuxData/MonfCrops.csv", stringsAsFactors = F) 
monf_cat <- unique(crops$Monf_Name)

fn <- Sys.glob("D:/Dp_world/CropSuit/byEcoID/*.tif")
suit <- rast(fn)

RIDs <- as.numeric(sub("D:/Dp_world/CropSuit/byEcoID/ID", "",
                       sub(".tif", "", fn)))
dirn <- "D:/Dp_world/CropSuit/byMonfCat"
dir.create(dirn)

for(i in 1:length(monf_cat)){
  ridi <- crops$RID[crops$Monf_Name == monf_cat[i]]
  ri <- which(RIDs %in% ridi) 
  # max values for the crop category 
  app(suit[[ri]], max, na.rm = T,
      filename = file.path(dirn, paste0(monf_cat[i], ".tif")),
      overwrite = T, 
      wopt = list(names = monf_cat[i], filetype = "GTiff"))
}


