# methods for sampling presence (only) observations from abundance data
## Arguments
## d = a data.frame or data.table
## crop = name of column with relative crop abundance data
## method = method for sampling presence (integer 1:7)
## npr = total number of presence observations to be sampled. 
## totcl = name of column with total cropland per cell
## colk = name of column with folds (for cross validation)
## k = fold for which suitability will be estimated (fold k won't be sampled)

fpr <- function(d, crop, method, npr, totcl = "tot", colk = NULL, k = NULL){
  
  # probability of cells being identify as presence
  if(method == 1){
    prob <- d[[crop]] > 0 
  } else if(method == 2){
    prob <- d[[crop]]     
  } else if(method == 3){
    prob <- d[[crop]] * d[[paste0(crop, "_dqi")]]
  } else if(method == 4){
    prob <- d[[crop]] * d[[totcl]] / max(d[[totcl]])
  } else if(method == 5){         
    prob <- d[[crop]] * d[[totcl]]/max(d[[totcl]]) * d[[paste0(crop, "_dqi")]]
  } 
  
  # for cross validation
  if(!is.null(k)){
    prob <- ifelse(d[[colk]] == k, 0, prob) # set k fold to prob of 0  
  }
  # not enough presences?
  rpl = sum(prob > 0) < npr
  
  pr <- sample(1:nrow(d), npr, replace = rpl, prob = prob)
  
  return(pr)
}
