allocate <- function(d, abun_columns, suit_columns, start.alloc, as_rat){
  
  ncr <- length(abun_columns)
  crops <- gsub(".abun", "", abun_columns)

  # suitability normalized to sum 1 for each cell (row)
  s1 <- paste0(crops, ".s1")
  d[, ss:= Reduce(`+`, .SD), .SDcols = suit_columns]
  d[, (s1):= lapply(.SD, function(x) x/ss), .SDcols = suit_columns]
  d[, ss:= NULL]
  
  # cropland per cell
  d[, cl:= Reduce(`+`, .SD), .SDcols = abun_columns]
  
  # total crop area per crop
  crop_area <- colSums(d[, ..abun_columns])
  names(crop_area) <- crops
  
  # create columns to add allocated area to crop
  allo_columns <- gsub(".suit", ".allo", suit_columns)
  d[, (allo_columns):= 0]
  
  # remaining area
  d[, ra:= cl]
  
  # set loops
  maxiter <- 100
  todo <- crop_area
  continue <- TRUE
  
  # allocation loop ----------
  while(continue){
    
    for(iter in 1:maxiter) {
      # fraction of cropland area per cell to be allocated
      f <- max(iter/maxiter, start.alloc) 
      
      # randomize crop order
      crops <- sample(crops)
      
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
          d[, (ja):= d[[ja]] + toadd] 
        }
        
        # remaining cropland area per cell after crop j
        d[, ra:= cl - Reduce(`+`, .SD), .SDcols = allo_columns]
      }
      
      # what is left to do
      todo <- crop_area - colSums(d[, ..allo_columns])
    }
    # has it converged?
    continue <- abs(todo) < crop_area * 0.01
    
    # if not, remove some cropland and start over
    if(continue){
      for (j in crops) {
        if(todo[j] > 0) next        # don't remove from under-allocated crops
        js <- paste0(j, ".suit")
        js1 <- paste0(j, ".s1")
        ja <- paste0(j, ".allo")
        
        # in order of relative suitability
        setorderv(d, js1, 1L)
        
        # cell with some j crop area
        i <- d[, .I(.SD > 0), .SDcols = ja]
        
        # crop area per cell to remove: based on relative suitability ranking  
        torm <- d[[ja]][i] * (1 - (1:length(i))/length(i))
        
        # remove
        d[i, (ja):= .SD - torm, .SDcols = ja] 
      }
      
      # remaining cropland area per cell after removing some cropland
      d[, ra:= cl - Reduce(`+`, .SD), .SDcols = allo_columns]
      # update to do
      todo <- crop_area - colSums(d[, ..allo_columns])
    }
  }
  return(d)
}
