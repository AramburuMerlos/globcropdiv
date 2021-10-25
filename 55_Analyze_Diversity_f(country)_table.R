library(magrittr)
library(data.table)
library(flextable)
library(officer)
library(officedown)

if(system('hostname', TRUE) == "ESP-RH-9891"){
  setwd("D:/globcropdiv/")
} else if(system('hostname', TRUE) == "LAPTOP-ST129J47"){ 
  setwd("G:/My Drive/globcropdiv/")
} # else if { ... 

tdir <- "G:/My Drive/globcropdiv/Tables"
dir.create("tdir", F, T)


# import data
d <- fread("OutData/D_by_country.csv")

countries <- geodata::country_codes()
d <- merge(d, countries[, c("NAME", "ISO3", "continent")],
           by.x= "country", by.y = "ISO3")

# subset and reorder
d <- d[NAME != "Caspian Sea"]
#d <- d[tot.area > 20000, ]
setorderv(d, "tot.area", -1L)


# keep only alpha diversity
Dg_cols <- names(d)[names(d) %like% "Dg"]
d[, (Dg_cols):= NULL]
# remove potential diversity
pot_cols <- names(d)[names(d) %like% "pot"]
d[, (pot_cols):= NULL]

# attainable diversity average
d[, att_Da_avg:= rowMeans(cbind(att_Da_sdm, att_Da_eco), na.rm = T)]

# compute diversity gap
d[, gap_Da_avg:= (att_Da_avg - act_Da)/att_Da_avg * 100]


# area to text
d[tot.area >= 1e6, area:= paste(format(tot.area/1e6, digits = 2), "M")]
d[tot.area < 1e6 & tot.area >= 1e3, area:= paste(round(tot.area/1e3), "k")]
d[tot.area < 1e3, area:= as.character(round(tot.area))]

db <- d[tot.area > 1e5]


d[, continent:= gsub("North", "N.", continent)]
d[, continent:= gsub("South", "S.", continent)]

col.ord <- c("NAME", "continent", 
             "act_Da","att_Da_avg", "gap_Da_avg",
             "area")

setcolorder(d, col.ord)

db <- d[tot.area > 1e5]


d[, c("att_Da_eco", "att_Da_sdm", "country", "tot.area"):= NULL]
db[, c("att_Da_eco", "att_Da_sdm", "country", "tot.area"):= NULL]



# create objects
col.list <- c("Country", "Continent",
           "Current", "Attainable", "Gap (%)", 
           "Cropland (ha)")
names(col.list) <- col.ord

head1 <-  c("Country", "Continent",
            "Local Diversity Average", 
            "Cropland (ha)")

brdr1 <- fp_border(color = "#666666", width = 1.5)
brdr2 <- fp_border(color = "#666666", width = .75)

d[NAME == "Democratic Republic of the Congo", NAME:= "DRC"]




flextable(d) %>% 
  set_header_labels(values = col.list) %>% 
  add_header_row(colwidths = c(1, 1, 3, 1), values = head1) %>% 
  merge_v(j = c(1:2, 6), part = "header") %>% 
  theme_zebra(odd_header = "transparent", even_header = "transparent") %>%
  hline(i = 1, j= 3:5, border = brdr2, part = "header") %>% 
  hline(i = 2, border = brdr1, part = "header") %>%
  hline(i = 1, j= c(1:2,6), border = brdr1, part = "header") %>% 
#  compose(1, 3, as_paragraph(as_sup("c"), as_i("D")), "header") %>% 
#  compose(1, 4, as_paragraph(as_sup("a"), as_i("D")), "header") %>% 
  align(align = "center") %>% 
  align(j = c(1,2), align = "left") %>% 
  align(align = "center", part = "header") %>% 
  align(j = c(1,2), align = "left", part = "header") %>% 
  colformat_double(j = c(5), digits = 0) %>% 
  colformat_double(j = c(3:4), digits = 1) %>% 
  vline(j = c(1,3:4), border = brdr2) %>% 
  vline(j = c(2,5), border = brdr2) %>% 
  border_outer(brdr1) %>% 
  fontsize(size = 10, part = "body") %>% 
  fontsize(size = 11 , part = "header") %>%  
  set_table_properties(layout = "autofit") %>% 
  autofit(add_w = 0.08, add_h = 0.08) %>% 
  set_caption("Sup. Table 1. Country local diversity averages and gaps. Local average refers to the weighted average of all 5 arc minutes resolution cropland cells inside a country. Attainable diversity is the average attainable diversity estimated with Ecocrop and Spatial Distribution Models. Diversity gap is the difference between attainable and current diversities, expressed as percentage of the attainable diversity. In cropland area, M refers to millions ha and k to thousands of ha.", style = "Normal") %>%  
              # Countries with less than 20,000 ha were ignored.", style = "Normal") %>% 
  save_as_docx(path = file.path(tdir, "Country_D.docx"),  
               pr_section = prop_section(page_margins = page_mar(.5,.5,.5,.5),
                                         page_size = page_size(8.5,11)))

