library(terra)

nwd <- nchar(getwd())
if(substr(getwd(),  nwd - 10, nwd) != "globcropdiv") stop("See 0000_wd.R")

path <- "D:/globcropdiv"
path <- ifelse(dir.exists(path), path, "OutData") 

# for crop species with multiples entries in the EcoCrop Data base, 
# keep the highest suitability value for the species. 
crops <- read.csv("AuxData/MonfCrops.csv", stringsAsFactors = F) 
monf_cat <- unique(crops$Monf_Name)

fn <- Sys.glob(file.path(path, "EcocropSuit/byEcoID/*.tif"))
suit <- rast(fn)

RIDs <- as.numeric(sub(file.path(path, "EcocropSuit/byEcoID/ID"), "",
                       sub(".tif", "", fn)))
dirn <- file.path(path, "EcocropSuit/byMonfCat")
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


