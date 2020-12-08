library(terra)
library(data.table)
library(corrplot)

nwd <- nchar(getwd())
if(substr(getwd(),  nwd - 10, nwd) != "globcropdiv") warning("See 0000_wd.R")

if(dir.exists("D:/globcropdiv")){
  wcpath <- "D:/WorldClim/2.1/wc5min"
  sgpath <- "D:/SoilGrids"
  gdpath <- "D:/globcropdiv"
  aqpath <- "D:/AQUASTAT"
  cpath <- "D:/WorldCropAbundance"
} else {
  wcpath <- "InData/WorldClim/2.1/wc5min"
  sgpath <- "InData/SoilGrids"
  aqpath <- "InData/AQUASTAT"
  gdpath <- "OutData"
  cpath <- "InData/WorldCropAbundance"
}

# Load Predictors and Mask Cropland --------------
preds <- c(rast(Sys.glob(file.path(wcpath, "bioc/*.tif"))), # WorldClim Bioclim
           rast(file.path(wcpath, "extra/GDD.tif")),
           rast(Sys.glob(file.path(wcpath, "extra/AI?q.tif"))),
           rast(file.path(wcpath, "extra/AI.tif")),  
           rast(file.path(wcpath, "extra/PET_seasonality.tif")),  
           rast(file.path(sgpath, "phh2o/phh2o_0-15cm_mean_5min.tif")), 
           rast(file.path(aqpath, "gmia_v5_aei_pct.asc"))  # irrigation
)

cmask <- rast(file.path(cpath, "Total_Cropland_ha.tif"))
preds <- mask(preds, cmask)

# Get Values and rename preds
DTp <- as.data.table(values(preds))
DTp <- DTp[complete.cases(DTp),]
oldn <- names(DTp)
setnames(DTp, 
         oldn[grepl("wc2.1", oldn)], 
         gsub("wc2.1_5m_bio_", "bc", oldn[grepl("wc2.1", oldn)]))
setnames(DTp, 
         c("Aridity_Index", "PET_seasonality", "pH0_15", "gmia_v5_aei_pct"),
         c("AI", "PETs", "pH", "Irr"))

# Correlation Matrix to identify highly correlated variables
M <- cor(DTp)

corrplot.mixed(M, order = "AOE")

# PCA to see which variables explain most of the first dimensions and to keep enough variables to explain most of the variation. Only climatic variables
d <- DTp[,-c("Irr", "pH")]
# PCA analysis
pca <- prcomp (d, scale = TRUE)
summary(pca)

# first 2 PC
plot(t(cor(pca$x[,1:2], d)), type = "n",
     ylim = c(-1,1), xlim = c(-1,1),
     main = "The variables in PC1 and PC2",
     pch = 20)
abline (h=0,lty=3)
abline (v=0,lty=3)
arrows(0,0,t(cor(pca$x[,1], d)),t(cor(pca$x[,2], d)))
text(t(cor(pca$x[,1], d)), t(cor(pca$x[,2], d)), names(d), pos=3 , cex = 0.8) 

# 3rd and 4th
plot(t(cor(pca$x[,3:4], d)), type = "n",
     ylim = c(-1,1), xlim = c(-1,1),
     main = "The variables in PC3 and PC4",
     pch = 20)
abline (h=0,lty=3)
abline (v=0,lty=3)
arrows(0,0,t(cor(pca$x[,3], d)),t(cor(pca$x[,4], d)))
text(t(cor(pca$x[,3], d)), t(cor(pca$x[,4], d)), names(d), pos=3 , cex = 0.8) 

# removed predictors
rmp <- c("bc1", "bc11", "bc6", "bc10", "bc4", "bc13", "bc14", "AIcq", "AIhq")
DTrmp <- DTp[, ..rmp]

# selected predictors
DTp[, (rmp):= NULL]

# min porcentage of variance explained of removed predictors with selected predictors
min(sapply(DTrmp, function(x) summary(lm(x ~ ., data = DTp))$r.squared))

# create file with selected file names, path name and short name
snpred <- names(DTp)
npred <- gsub("bc", "wc2.1_5m_bio_", snpred)
fn <- c("PET_seasonality", "phh2o_0-15cm_mean_5min", "gmia_v5_aei_pct")
sn <- c("PETs", "pH", "Irr")
for(i in 1:length(fn)) npred <- gsub(sn[i], fn[i], npred)

outd <- data.table("short_name" = snpred, "name" = npred)
outd[, "path" := ifelse(grepl("^wc", name), "WorldClim/2.1/wc5min/bioc",
                        ifelse(grepl("^phh2o", name), "SoilGrids/phh2o",
                               ifelse(grepl("^gmia", name), "AQUASTAT",
                                      "WorldClim/2.1/wc5min/extra")))]
outd[, "file_name":= paste0(file.path(path, name), 
                            ifelse(grepl("^gmia", name), ".asc", ".tif"))]

fwrite(outd, "AuxData/SelectedPreds.csv")