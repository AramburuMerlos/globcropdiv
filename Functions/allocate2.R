allocate <- function(d, abun_columns, suit_columns, start.alloc, as_rat){
  ncr <- length(abun_columns)
  crops <- gsub(".abun", "", abun_columns)

  # suitability normalized to 1 per cell (row)
  s1 <- paste0(crops, ".s1")
  d[, ss:= Reduce(`+`, .SD), .SDcols = suit_columns]
  d[, (s1):= lapply(.SD, function(x) x/ss), .SDcols = suit_columns]
  
  
  # cropland per cell 
  d[, cl:= Reduce(`+`, .SD), .SDcols = abun_columns]
  
  # total crop area per crop 
  crop_area <- colSums(d[, ..abun_columns])
  names(crop_area) <- crops
  
  crops <- d[,..s1] * d$cl %>%
    colSums() %>%
    `-`(crop_area, .) %>%
    order(decreasing = T) %>%
    crops[.]
  
  # create columns to add crop allocated area and cropland proportion
  area_columns <- gsub(".suit", ".area", suit_columns)
  prop_columns <- gsub(".suit", ".prop", suit_columns)
  d[, (area_columns):= 0]
  d[, (prop_columns):= 0]
  
  # remaining area
  d[, ra:= cl]
  
  # initial allocation ------------
  for(j in crops){
    maxiter <- 5
    iter <- 0
    
    while(iter < maxiter) {
      iter = iter + 1
      
      js <- paste0(j, ".suit")
      js1 <- paste0(j, ".s1")
      ja <- paste0(j, ".area")
      jp <- paste0(j, ".prop")
      
      if(iter == 1){
        # factor
        f <- crop_area[j] / sum(d[[js1]] * d$ra)
        # allocated crop proportion
        d[, (jp):= pmin(d[[js1]] * f, 1)]
        # allocated crop area
        d[, (ja):= d[[jp]] * ra]
        
      } else {
        rem <- crop_area[j] - sum(d[[ja]])
        i = which(d[[js1]] < 1)
        
      }
        # update all relative suitabilities to honor jp
        d[, (js1):= d[[jp]]]
        other <- s1[s1!=js1]
        d[, ss:= Reduce(`+`, .SD), .SDcols = other]
        for(k in other) d[, (k):= d[[k]] * (1-d[[js1]])/d$ss]
        
     
      
      
      # update remaining area
      d[, ra:= cl - Reduce(`+`, .SD), .SDcols = area_columns]
      
      # done with this crop?
      if(sum(d[[ja]]) >= crop_area[j]) break
    }
  }  
  
  
  
    
    # fraction of cropland area per cell to be allocated
    f <- max(iter/maxiter, start.alloc) 
    
    crops <- sample(crops)  # randomize order (maybe narrow niche first?)
    
    for (j in crops) {
      if(todo[j] == 0) next
      js <- paste0(j, ".suit")
      js1 <- paste0(j, ".s1")
      ja <- paste0(j, ".allo")
      
      # in order of suitability
      setorderv(d, js, -1L)
      
      # area per cell that could be allocated to j in this iter (avoid 0 suit)
      toadd <- pmin(f * d$ra, (d$cl * d[[js1]] * as_rat) - d[[ja]])
      
      # is this the last round for j?
      if(sum(toadd) >= todo[j]){
        i.stop <- Position(function(x) x >= todo[j], cumsum(toadd))
        d[1:(i.stop - 1), (ja):= .SD + toadd[1:(i.stop - 1)], .SDcols = ja]
        d[i.stop, (ja):= d[[ja]][i.stop] + (crop_area[j] - sum(d[[ja]]))]
      } else {
        d[, (aj):= d[[aj]] + toadd] 
      }
      
      # remaining cropland area per cell after crop j
      d[, ra:= cl - Reduce(`+`, .SD), .SDcols = allo_columns]
    }
    
    # what is left todo
    todo <- crop_area - colSums(d[, ..allo_columns])
    
    # are we done?
    if (all(todo == 0)) break
  }
  
  
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
  
  
  
  
return(d)
}

