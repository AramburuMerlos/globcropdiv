library(magrittr)
library(data.table)
library(flextable)
library(officer)
library(officedown)

if(system('hostname', TRUE) %in% c("ESP-RH-9891", "LAPTOP-ST129J47")){
  setwd("D:/globcropdiv/")
} # else if { ... 


tdir <- "G:/My Drive/globcropdiv/Tables"
dir.create(tdir, F, T)


# import data
d <- fread("OutData/D_by_country.csv")

countries <- geodata::country_codes()
d <- merge(d, countries[, c("NAME", "ISO3", "continent")],
           by.x= "country", by.y = "ISO3")

# subset and reorder
d <- d[NAME != "Caspian Sea"]
#d <- d[tot.area > 20000, ]
setorderv(d, "tot.area", -1L)


# remove potential diversity
pot_cols <- names(d)[names(d) %like% "pot"]
d[, (pot_cols):= NULL]

# attainable diversity average
d[, att_Da_avg:= rowMeans(cbind(att_Da_sdm, att_Da_eco), na.rm = T)]
d[, att_Dg_avg:= rowMeans(cbind(att_Dg_sdm, att_Dg_eco), na.rm = T)]


# compute diversity gap
d[, gap_Da_avg:= (att_Da_avg - act_Da)/att_Da_avg * 100]
d[, gap_Dg_avg:= (att_Dg_avg - act_Dg)/att_Dg_avg * 100]


# area to text
d[tot.area >= 1e6, area:= paste(format(tot.area/1e6, digits = 2), "M")]
d[tot.area < 1e6 & tot.area >= 1e3, area:= paste(round(tot.area/1e3), "k")]
d[tot.area < 1e3, area:= as.character(round(tot.area))]

db <- d[tot.area > 1e5]

# interesting values ------
setorder(d, act_Da)
d[,.(NAME, act_Da, area)]

setorder(db, act_Da)
db[1, .(NAME, act_Da)]
db[act_Da > 12, .(NAME, act_Da)]
db[act_Da > 4 & act_Da < 8, .N]
db[act_Da < 3, .(NAME, act_Da, area)]

db[act_Dg  > 20, .N]
db[act_Dg  > 12, .N]/db[, .N]
db[act_Da > 4 & act_Da < 8, .N]
db[act_Da < 3]

setorder(db, gap_Da_avg)
db[1:5, .(NAME, gap_Da_avg)]
db[gap_Da_avg < 50, .N]

db[gap_Da_avg > 70, .N]/db[, .N]

setorder(db, gap_Dg_avg)
db[1:10, .(NAME, gap_Dg_avg)]
db[gap_Dg_avg < 50, .N]
db[gap_Da_avg > 70, .N]/db[, .N]
hist(db$gap_Dg_avg)
hist(db$gap_Da_avg)

plot(d$gap_Dg_avg, d$gap_Da_avg, xlim = c(0,100), ylim = c(0,100))
abline(a = 0, b = 1, col = "red")
t.test(d$gap_Dg_avg, d$gap_Da_avg, alternative = "less", paired = TRUE)

# end ------

d[, continent:= gsub("North", "N.", continent)]
d[, continent:= gsub("South", "S.", continent)]

col.ord <- c("NAME", "continent", 
             "act_Da","att_Da_avg", "gap_Da_avg",
             "act_Dg","att_Dg_avg", "gap_Dg_avg",
             "area")

setcolorder(d, col.ord)


d[, c("att_Da_eco", "att_Da_sdm"):= NULL]
d[, c("att_Dg_eco", "att_Dg_sdm"):= NULL]
d[, c("country", "tot.area"):= NULL]

head(d)

# create objects
col.list <- c(
  "Country", "Continent", 
  rep(c("cD", "aD", "Dg(%)"), 2), 
  "Cropland (ha)"
)

names(col.list) <- col.ord

head1 <-  c("Country", "Continent",
            paste0("D", as.character("\u03B1")),
            paste0("D", as.character("\u03B3")),
            "Cropland (ha)")

brdr1 <- fp_border(color = "#666666", width = 1.5)
brdr2 <- fp_border(color = "#666666", width = .75)

d[NAME == "Democratic Republic of the Congo", NAME:= "DRC"]




flextable(d) %>% 
  set_header_labels(values = col.list) %>% 
  add_header_row(colwidths = c(1, 1, 3,3, 1), values = head1) %>% 
  merge_v(j = c(1:2, 9), part = "header") %>% 
  theme_zebra(odd_header = "transparent", even_header = "transparent") %>%
  hline(i = 1, j= 3:8, border = brdr2, part = "header") %>% 
  hline(i = 2, border = brdr1, part = "header") %>%
  hline(i = 1, j= c(1:2,9), border = brdr1, part = "header") %>% 
#  compose(1, 3, as_paragraph(as_sup("c"), as_i("D")), "header") %>% 
#  compose(1, 4, as_paragraph(as_sup("a"), as_i("D")), "header") %>% 
  align(align = "center") %>% 
  align(j = c(1,2), align = "left") %>% 
  align(align = "center", part = "header") %>% 
  align(j = c(1,2), align = "left", part = "header") %>% 
  colformat_double(j = c(5,8), digits = 0) %>% 
  colformat_double(j = c(3:4,6:7), digits = 1) %>% 
  vline(j = c(1,3:4,6:8), border = brdr2) %>% 
  vline(j = c(2,5,9), border = brdr2) %>% 
  border_outer(brdr1) %>% 
  fontsize(size = 10, part = "body") %>% 
  fontsize(size = 11 , part = "header") %>% 
  italic(j = c(3:8), part = "header") %>% 
  set_table_properties(layout = "autofit") %>% 
  autofit(add_w = 0.08, add_h = 0.08) %>% 
  set_caption("Sup. Table 1. Country local diversity average (Dα) and total diversity (Dγ) for current (cD) and attainable (aD) crop species diversity levels and their corresponding diversity gap (Dg, %). Local diversity average refers to the weighted average of all 86km2 cropland cells inside a country. aD is the average attainable diversity estimated with two methods. Dg is the difference between aD and cD, expressed as a percentage of the aD. For cropland area, M refers to millions ha and k to thousands of ha.", style = "Normal") %>%  
              # Countries with less than 20,000 ha were ignored.", style = "Normal") %>% 
  save_as_docx(path = file.path(tdir, "Country_D.docx"),  
               pr_section = prop_section(page_margins = page_mar(.5,.5,.5,.5),
                                         page_size = page_size(8.5,11)))

