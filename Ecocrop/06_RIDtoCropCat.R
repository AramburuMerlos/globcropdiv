library(magrittr)
library(terra)

if(!grepl("globcropdiv$", getwd())) warning("See 0000_wd.R")

# for crop species with multiples entries in the EcoCrop Data base, 
# keep the highest suitability value for the species. 
d <- read.csv("AuxData/EcocropSpeciesWithAbundanceData.csv", 
              stringsAsFactors = F) 
crops <- d%$% ifelse(Monf_Name == "", SPAM_Name, Monf_Name) %>% unique

suit <- "OutData/EcocropSuit/byEcoID/*.tif" %>% Sys.glob %>% rast

RIDs <- names(suit) %>% sub("ID", "",.) %>% as.numeric

outdir <- file.path("OutData/EcocropSuit/byCrop")
dir.create(outdir, F, T)

for(i in 1:length(crops)){
  ridi <- d$RID[d$Monf_Name == crops[i] | d$SPAM_Name == crops[i]]
  ri <- which(RIDs %in% ridi) 
  # max values for the crop category 
  app(suit[[ri]], max, na.rm = T,
      filename = file.path(outdir, paste0(crops[i], ".tif")),
      overwrite = T, 
      wopt = list(names = crops[i], filetype = "GTiff"))
}


