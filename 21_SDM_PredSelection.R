library(magrittr)
library(terra)
library(data.table)
library(corrplot)

if(!grepl("globcropdiv$", getwd())){
  if (system('hostname', TRUE) %in% c("ESP-RH-9891", "LAPTOP-ST129J47")) { 
    setwd("G:/My Drive/globcropdiv/")
  } # else if { ... 
}


# Load Data -----------
# Prepare data table
DT <- "InData/TotalCropland.tif" %>% 
  rast() %>% 
  values() %>% 
  is.na() %>% 
  `!`() %>% 
  which() %>% 
  data.table(cell = .)
  
# Predictors
wcdir <- "InData/WorldClim/2.1/wc5min"

r_preds <- c(rast(Sys.glob(file.path(wcdir, "bioc/*.tif"))), # WorldClim Bioclim
             rast(file.path(wcdir, "extra/GDD.tif")),
             rast(Sys.glob(file.path(wcdir, "extra/AI?q.tif"))),
             rast(file.path(wcdir, "extra/AI.tif")),  
             rast(file.path(wcdir, "extra/PET_seasonality.tif")),  
             rast("InData/SoilGrids/phh2o/phh2o_0-15cm_mean_5min.tif"), 
             rast("InData/AQUASTAT/gmia_v5_aei_pct.asc")  # irrigation
)

# extract values and rename columns
oldn <- names(r_preds)
DT[, (oldn):= extract(r_preds, cell)]
DT[, cell:= NULL]
setnames(DT, 
         oldn[grepl("wc2.1", oldn)], 
         gsub("wc2.1_5m_bio_", "bc", oldn[grepl("wc2.1", oldn)]))
setnames(DT, 
         c("Aridity_Index", "PET_seasonality", "pH0-15HR", "gmia_v5_aei_pct"),
         c("AI", "PETs", "pH", "Irr"))

# Correlation Matrix to identify highly correlated variables
M <- cor(DT)

corrplot.mixed(M, order = "AOE")

# PCA to see which variables explain most of the first dimensions and to keep enough variables to explain most of the variation. Only climatic variables
d <- DT[,-c("Irr", "pH")]
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
DTrmp <- DT[, ..rmp]

# selected predictors
DT[, (rmp):= NULL]

# percentage of removed predictors variance explained by selected predictors
sapply(DTrmp, function(x) summary(lm(x ~ ., data = DT))$r.squared)

# create file with selected file names, path name and short name
snpred <- names(DT)
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