# install.packages(c("meteor", "terra"))
# remotes::install_github("cropmodels/Recocrop")
library(Recocrop)
library(data.table)
dirname <- "AuxData"
dir.create(dirname)

cn <- ecocropPars() # crop names
setDT(cn)
cn[, RID:= 1:.N]
setcolorder(cn, "RID")
setorder(cn, "SCIENTNAME")
cn
fwrite(cn, file.path(dirname, "Ecocrop_names.csv"))

# dowload crop data from FAOSTAT (world totals)
# link: http://www.fao.org/faostat/en/#data/QC (retrieved 08/10/2020)
fao <- fread(file.path(dirname, "FAOSTAT_world_data.csv"))
setnames(fao, names(fao), gsub(" ", "_", names(fao)))
fcrops <- unique(fao[,.(Item_Code, Item)])
setorder(fcrops, Item_Code)
fwrite(fcrops, file.path(dirname, "FAOcrops.csv"))

#############################

# match FAO crops with ECOCROP crop species externally.
# only crops included in FAOSTAT will be considered. Species level
# see Ecicrios_metadata.txt for a description of the following table
ecrops <- fread(file.path(dirname, "Ecocrops.csv"))

# check if Row ID and names match
setorder(cn, RID)
setorder(ecrops, RID)
mapply(function(x,y) all.equal(x,y), 
       cn,
       ecrops[,.(RID, NAME, SCIENTNAME)])

# which FAO crops aren't represented in the selected crops
fcrops[!fcrops$Item_Code %in% ecrops$FAO_Code &
         !fcrops$Item_Code %in% ecrops$FAO_Code2 &
         !fcrops$Item_Code %in% ecrops$FAO_Code3]
# Triticale (97) isn't in Ecocrop database
# Grain, mixed (103) doesn't specify the species included

logcols <- c("FAOSTAT", "crop", "food", "fiber", "fodder", "forestry", "industry", "plant")
for(j in logcols) set(ecrops, j = j, value = as.logical(ecrops[[j]]))
# selected crop species: cultivated plants for food included in FAOSTAT 
sc <- ecrops[FAOSTAT & crop & food & plant,]
nrow(sc)

head(sc)
fwrite(sc, file.path(dirname, "selsp.csv"))
