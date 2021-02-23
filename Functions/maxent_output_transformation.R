# return transformed cloglog values after calibrating c to match mean abundance
# (after Guillera-Arroita et al., 2014)
clogtr <- function(mean_abun, lp, maxiter = 500){
  
  # cloglog transformation function
  .fcloglog <- function(c,lp) 1 - exp(-exp(lp + c))
  
  # Steven Phillips starting values
  cL <- ifelse(mean(.fcloglog(c = -5, lp = lp)) > mean_abun, -20, -5)
  cH <- ifelse(mean(.fcloglog(c = 5, lp = lp)) < mean_abun, 20, 5)
  
  # search the value of c that makes mean(cloglog) == mean(abun)
  # and return transformed values
  ii <- 0
  continue <- TRUE
  while(continue) {
    cM <- (cH+cL)/2   
    new_est <- .fcloglog(cM, lp)
    dif <- mean(new_est) - mean_abun
    continue <- abs(dif) > 0.01
    if(continue){
      if(dif > 0) cH <- cM else cL <- cM 
      ii <- ii + 1
      if(ii == maxiter){
        continue <- FALSE
        warning("c value didn't converge, increase maxiter and/or review lp")
      } 
    }
  }
return(new_est)
}
  




