abun_columns <- ca
suit_columns <- cs
rm(ca, cs, j)

# BACKUP AND DEBUG --------  
#db <- data.table::copy(d)
d <- data.table::copy(db)

uj = uac[1]
oj = oac[1]


rm(list = ls()[ls()!= "db"])


###### CUT from allocate Function ########

# fill remaining cropland -----------
# with most suitable crop 

ia <- d[,.I[ra > 0]]
for(j in crops) {
  js1 <- paste0(j, ".s1")
  ja <- paste0(j, ".allo")
  d[ia, (ja):= d[[ja]][ia] + d[[js1]][ia] * d$ra[ia]]
}

d[, ra:= cl - Reduce(`+`, .SD), .SDcols = allo_columns]

# reallocate cropland ---------------------
# to match total crop areas

todo <- crop_area - colSums(d[, ..allo_columns])

# under allocated crops
uac <- todo[todo > 0]  %>% names()
uac <- uac[order(todo[uac], decreasing = T)]

# over allocated crops
oac <- todo[todo < 0] %>% names()
oac <- oac[order(todo[oac])]

# increase area of under allocated crops at the expense of over allocated
for(uj in uac){
  ujs <- paste0(uj, ".suit")
  uja <- paste0(uj, ".allo")
  
  for(oj in oac){
    # have we already used this extra crop area?
    if(todo[oj] >= 0) next
    
    ojs <- paste0(oj, ".suit")
    oja <- paste0(oj, ".allo")
    
    # order by best cells to do reallocation (based on suitability)
    d[, suit.rel:= d[[ujs]]-d[[ojs]]]
    setorderv(d, "suit.rel", -1L)
    d[, suit.rel:= NULL]
    
    # total area to be reallocated  
    tara <- min(todo[uj], -todo[oj])  
    
    # cells where area could be reallocated (suitable for uj, with oj area)
    i <- which(d[[ujs]] > 0 & d[[oja]] > 0)
    
    # none suitable cells for uj with crop oj?
    if(length(i) == 0) next
    
    # if there is enough oj area to reach tara
    if(sum(d[[oja]][i]) > tara){
      i.stop <- Position(function(x) x >= tara, cumsum(d[[oja]][i]))
      rest <- tara - sum(d[[oja]][i[1:(i.stop-1)]])
      
      d[i[1:(i.stop - 1)], (uja):= Reduce(`+`, .SD), .SDcols = c(uja, oja)]
      d[i[1:(i.stop - 1)], (oja):= 0]
      
      d[i[i.stop], (uja):= .SD + rest, .SDcols = uja]
      d[i[i.stop], (oja):= .SD - rest, .SDcols = oja]
      
    } else {  # if not enough, reallocate everything suitable with oj area
      d[i, (uja):= Reduce(`+`, .SD), .SDcols = c(uja, oja)]
      d[i, (oja):= 0]
    }
    
    # update to do for crop uj
    todo[uj] <- crop_area[uj] - sum(d[[uja]])
    # update to do for crop oj
    todo[oj] <- crop_area[oj] - sum(d[[oja]])
  }
}



## step 2 ----
# use AN intermediary crop to solve unbalances

# under allocated crops
uac <- todo[todo > 0]  %>% names()
uac <- uac[order(todo[uac], decreasing = T)]

# over allocated crops
oac <- todo[todo < 0] %>% names()
oac <- oac[order(todo[oac])]

# possible intermediary crops
iac <- todo[todo == 0] %>% names()
iac.suit <- paste0(iac, ".suit")
iac.allo <- paste0(iac, ".allo")

# transfer area of over allocated crops to the intermediary
for(uj in uac){
  ujs <- paste0(uj, ".suit")
  uja <- paste0(uj, ".allo")
  
  for(oj in oac){
    # have we already used this extra crop area?
    if(todo[oj] >= 0) next
    
    ojs <- paste0(oj, ".suit")
    oja <- paste0(oj, ".allo")
    
    # choose intermediary
    sapply(d[,..iac.suit], function(x) cor(x,d[[ujs]]))
    sapply(d[,..iac.allo], function(x) cor(x,d[[oja]]))
    
    # order by best cells to do reallocation (based on suitability)
    d[, suit.rel:= d[[ijs]]-d[[ojs]]]
    setorderv(d, "suit.rel", -1L)
    d[, suit.rel:= NULL]
    
    # cells where area can be reallocated (suitable for maize, with oj area)
    i <- which(d[[ijs]] > 0 & d[[oja]] > 0)
    
    # none suitable cells for uj with crop oj?
    if(length(i) == 0) next
    
    if(sum(d[[oja]][i]) > -todo[oj]){
      i.stop <- Position(function(x) x >= tara, cumsum(d[[oja]][i]))
      rest <- tara - sum(d[[oja]][i[1:(i.stop-1)]])
      
      d[i[1:(i.stop - 1)], (uja):= Reduce(`+`, .SD), .SDcols = c(uja, oja)]
      d[i[1:(i.stop - 1)], (oja):= 0]
      
      d[i[i.stop], (uja):= .SD + rest, .SDcols = uja]
      d[i[i.stop], (oja):= .SD - rest, .SDcols = oja]
      
    } else {  # if not enough, reallocate everything suitable with oj area
      d[i, (uja):= Reduce(`+`, .SD), .SDcols = c(uja, oja)]
      d[i, (oja):= 0]
    }
    
    # update to do for crop uj
    todo[uj] <- crop_area[uj] - sum(d[[uja]])
    # update to do for crop oj
    todo[oj] <- crop_area[oj] - sum(d[[oja]])
  }
  return(d)
}

}