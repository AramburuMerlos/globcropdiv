# d = a data.table with crop suitability and total cropland per cell (wide)
# crops (character) = a vector with crop names
# crop_area (named numeric) = a vector with crop areas to be allocated
# scols (character) = a vector with column names with suitability data
# tccol (character) = name of column with total cropland area per cell 

# Crops, crop_area and scols order must match

allocate <- function(d, crops, crop_area,  
                     scols = paste0("s.", crops),
                     tccol = "tcl",
                     showProgress = FALSE){
  
  if(!identical(crops, names(crop_area))) {
     stop ("crops and crop_area names must match") 
  }
  if(!all(mapply(grepl, crops, scols))) {
    stop ("crops and scols names must match") 
  }
  
  if(!"data.table" %in% class(d)) setDT(d)
  
  # set column names to default
  orig_scols <- scols
  scols <- paste0("s.", crops)
  setnames(d, orig_scols, scols)
  setnames(d, tccol, "tcl")
  
  if(round(sum(crop_area)/d[, sum(tcl)], 3) != 1){
    stop("cropland area is not equal to sum(crop_area) +/- 0.1%")
  }  
  
  # suitability normalized to sum 1 for each cell (row)
  ncols <- paste0("n.", crops)
  d[, ss:= Reduce(`+`, .SD), .SDcols = scols]
  d[, (ncols):= lapply(.SD, function(x) x/ss), .SDcols = scols]
  d[, ss:= NULL]
  
  # normalized suitability times cell cropland per crop
  nsxcl <- d[, lapply(.SD, function(x) sum(x*tcl, na.rm = T)), .SDcols = ncols]
  # crops relative abundance (factors)
  fact <- crop_area/unlist(nsxcl)
  #fact <- crop_area/mean(crop_area)
  
  # names of columns with allocated area
  acols <- paste0("a.", crops)

  # Initial allocation ---------
  for (j in crops) {
    nj <- paste0("n.", j)
    aj <- paste0("a.", j)
    # allocate area to j ~ relative suitability and crop relative abundance
    d[, (aj):= .SD * fact[j] * tcl, .SDcols = nj]
  }

  # row sums
  d[, rs:= Reduce(`+`, .SD), .SDcols = acols]
  
  # Remove extra cropland per cell
  d[rs > tcl, (acols):= lapply(.SD, function(x) x * tcl/rs), .SDcols = acols]
  
  # Remove extra crop_area per crop
  aar <- colSums(d[, ..acols])/crop_area # allocation area ratio
  for(aj in acols){
    if(aar[aj] > 1){
      d[, (aj):= .SD/aar[aj], .SDcols = aj]
    }
  }
  rm(aar)
  d[, rs:= NULL]
  
  # what is left to do
  todo <- crop_area - colSums(d[, ..acols])
  # remaining cropland per cell
  d[, ra:= tcl - Reduce(`+`, .SD), .SDcols = acols]
  d[ra < 0, ra:= 0] # to avoid rounding problems
  
  if(showProgress){
    mess <- paste0(round(sum(todo)/sum(crop_area) * 100, 2), 
                   "% area left to allocate")
    cat(paste(rep("\b", nchar(mess)), collapse = ""))
    cat(mess)
  }
  
  # allocation loop ----------
  maxiter <- 100
  continue <- TRUE
  
  while(continue){
    sum_todo_b4 <- sum(todo)
    for(iter in 1:maxiter) {
      # order crops by todo
      crops <- crops[order(todo, decreasing = TRUE)]
      
      for (j in crops) {
        if(todo[j] == 0) next
        sj <- paste0("s.", j)
        nj <- paste0("n.", j)
        aj <- paste0("a.", j)
        
        # fraction of cropland area per cell to be allocated. 
        # widespread crops: start with higher fraction
        f <- max(iter/maxiter, todo[j]/sum(todo)) 
        
        # in order of suitability
        setorderv(d, sj, -1L)
        
        # area per cell that could be allocated to j in this iter (avoid 0 suit)
        toadd <- pmin(f * d$ra, (d$tcl * d[[nj]] * 1000) - d[[aj]])
        toadd[toadd < 0] <- 0 # to avoid rounding problems
        
        # is this the last round for j?
        if(sum(toadd) >= todo[j]){
          i.stop <- min(which(cumsum(toadd) >= todo[j]))
          d[1:(i.stop - 1), (aj):= .SD + toadd[1:(i.stop - 1)], .SDcols = aj]
          d[i.stop, (aj):= d[[aj]][i.stop] + (crop_area[j] - sum(d[[aj]]))]
        } else {
          d[, (aj):= d[[aj]] + toadd] 
        }
        
        # remaining cropland area per cell after crop j
        d[, ra:= tcl - Reduce(`+`, .SD), .SDcols = acols]
        d[ra < 0, ra:= 0] # to avoid rounding problems
      }
      
      # what is left to do
      todo <- crop_area - colSums(d[, ..acols])
      if(showProgress){
        mess <- paste0(round(sum(todo)/sum(crop_area) * 100, 2), 
                       "% area left to allocate")
        cat(paste(rep("\b", nchar(mess)), collapse = ""))
        cat(mess)
      }
    }
    
    # has it converged?
    continue <- any(todo > (crop_area * 0.01))
    
    # if not, reallocate cropland to under allocated crops and start over
    if(continue){
      # under allocated crops in order of remaining area
      under <- names(which(todo[order(todo, decreasing = TRUE)] > 0))
      # fully allocated crops
      full <- names(todo[todo == 0])
      afull <- paste0("a.", full)
      for (uj in under) {
        suj <- paste0("s.", uj)
        auj <- paste0("a.", uj)
        setorderv(d, suj, -1L)
        # max area that could be reallocated (proportionally to suj)
        tora <- d[, ..afull][, .SD * d[[suj]]] 
        # needed area / max area ratio
        rat <- unname(todo[uj] / sum(tora))
        
        if(rat > 1) x = 1   # this is extremely unlikely, but just in case
        while(rat > 1){     
          x <- x + 1
          tora <- d[, ..afull][, .SD * (d[[suj]] ^ (1/x))] 
          rat <- todo[uj] / sum(tora)
          if(x == 10) rat <- 1
        }
        # reallocate
        for(afj in afull){
          set(d, j = afj, value = d[[afj]] - (tora[[afj]] * rat)) # remove
          set(d, j = auj, value = d[[auj]] + (tora[[afj]] * rat)) #add
        }
      }
      # update to do (it shouldn't change, though)
      todo <- crop_area - colSums(d[, ..acols])
      # update remaining area
      d[, ra:= tcl - Reduce(`+`, .SD), .SDcols = acols]
      d[ra < 0, ra:= 0] # to avoid rounding problems
    }
    
    if(sum(todo) >= sum_todo_b4){
      warning("no improvement")
      break
    }  
  }
  
  if(showProgress){
    cat("\n")  
  }
  
  
  d[, (ncols):= NULL]
  d[, ra:= NULL]
  
  # reallocate ------------
  ce0 <- .ce(d[,..acols], d[,..scols])
  
  if(showProgress){
    message("Allocation finished. \nCross Entropy = ", round(ce0))  
  }
  
  
  # all crop combinations pairs (without replacement, order doesn't matter)
  cc <- .comb(crops)
  setnames(cc, c("V1", "V2"), c("A", "B"))
  
  nfail <- 0
  ce_red <- 0
  
  if(showProgress){
    nsucc <- 0
  }
  
  while(nfail < nrow(cc)){
    cc <- cc[sample(1:nrow(cc)),]  # shuffle order
    # find pair of crops with room for reallocation
    for(i in 1:nrow(cc)){
      A <- cc[i,A]
      B <- cc[i,B]
      sA <- paste0("s.", A)
      sB <- paste0("s.", B)
      aA <- paste0("a.", A)
      aB <- paste0("a.", B)
      # expected area ratio between crops A and B for each cell
      d[, AB:= (d[[sA]] * fact[A]) / (d[[sB]] * fact[B])]
      # cells where area ratio is greater than expected
      iA2B <- which(d[[aA]]/d[[aB]] > d$AB)
      # cells where area ratio is lower than expected
      iB2A <- which(d[[aA]]/d[[aB]] < d$AB)
      if(length(iA2B) > 0 & length(iB2A) > 0) break
    }  
    
    # nothing more to reallocate?
    if(length(iA2B) == 0 | length(iB2A) == 0) break 
    
    # after some algebra, it turns out that:
    A2B <- (1/d$AB[iA2B] * d[[aA]][iA2B] - d[[aB]][iA2B])/(1 + 1/d$AB[iA2B])
    B2A <- (d$AB[iB2A] * d[[aB]][iB2A] - d[[aA]][iB2A])/(1 + d$AB[iB2A])
    
    # if sA == 0, then AB == 0, thus, reallocate all A area to B
    A2B[d$AB[iA2B] == 0] <- d[[aA]][iA2B[d$AB[iA2B] == 0]]
    # if sB == 0, then AB == Inf, thus, reallocate all B area to A
    B2A[is.infinite(d$AB[iB2A])] <- d[[aB]][iB2A[is.infinite(d$AB[iB2A])]]
    
    tA2B <- sum(A2B)
    tB2A <- sum(B2A)
    
    # if the total a2b is less than b2a, limit b2a; else, limit a2b
    if(tA2B < tB2A){
      i.stop <- min(which(cumsum(B2A) >= tA2B))
      B2A <- B2A[1:(i.stop - 1)]
      B2A <- if(i.stop == 1) tA2B else c(B2A, tA2B - sum(B2A))
      iB2A <- iB2A[1:i.stop]
    } else {
      i.stop <- min(which(cumsum(A2B) >= tB2A))
      A2B <- A2B[1:(i.stop - 1)]
      A2B <- if(i.stop == 1) tB2A else c(A2B, tB2A - sum(A2B))
      iA2B <- iA2B[1:i.stop]
    }
    
    # new area values for a and b
    naA <- c(d[[aA]][iA2B] - A2B, d[[aA]][iB2A] + B2A)     
    naB <- c(d[[aB]][iA2B] + A2B, d[[aB]][iB2A] - B2A)  
    iall <- c(iA2B, iB2A)
    
    # calculate Cross Entropy without reallocation for those cells
    ce_b4 <- .ce(d[iall, .SD, .SDcols = c(aA, aB)], 
                 d[iall, .SD, .SDcols = c(sA, sB)])
    # calculate Cross Entropy after reallocation for those cells
    ce_af <- .ce(data.table(naA, naB), 
                 d[iall, .SD, .SDcols = c(sA, sB)])           

    if(ce_af < ce_b4){
      d[iall, (aA):= naA]
      d[iall, (aB):= naB]
      ce_red <- ce_red + (ce_b4 - ce_af)
    }
    
    if(ce_af < (ce_b4 - ce0/1e4)){
      nfail <- 0
      if(showProgress) nsucc <- nsucc + 1
    } else {
      nfail <- nfail + 1
    }
    if(showProgress){
      if(nsucc == 2e5){
        prog <- round((1 - .ce(d[,..acols], d[,..scols])/ce0) * 100, 2)
        mess <- paste0("Cross Entropy Reduced in ", prog, "%")
        cat(paste(rep("\b", nchar(mess)), collapse = ""))
        cat(mess)
        nsucc <- 0
      }
    }
  }
  
  if(showProgress){
    cat("\n")  
  }
  
  d[, AB:= NULL]
  
  # set column names back
  setnames(d, scols, orig_scols)
  setnames(d, "tcl", tccol)
  
  return(d)
}




# Cross Entropy function
# x = data.frame with allocated area values
# y = data.frame with prior area values (prior = suitability)
.ce <- function(x, y, norm_by_row = TRUE, norm_by_col = FALSE){
  
  if(norm_by_row){
    x <- x/rowSums(x) # allocated area to cell fraction
    y <- y/rowSums(y) # area priors to cell fraction
  }
  
  if(norm_by_col){
    x <- x/colSums(x) # allocated area to allocated share
    y <- y/colSums(y) # area priors to share priors
  }
  
  x <- unlist(x, use.names = FALSE)
  y <- unlist(y, use.names = FALSE)
  
  if(length(x) != length(y)){ 
    stop("unequal vector lengths")
  }
  
  z <- x == 0 | y == 0
  x <- x[!z]
  y <- y[!z]
  
  sum(x * log(x)) - sum(x * log(y))
}




# Crop combinations (without replacement) function
.comb <- function(x){
  n <- length(x)
  l <- lapply(1:(n-1), function(i) data.table(rep(x[i], n-i), x[(i+1):n]))
  rbindlist(l)
} 




