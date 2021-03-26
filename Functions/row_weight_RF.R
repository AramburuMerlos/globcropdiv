# row weight in Random Forest. Based on total cropland and DQI

frw <- function(d, method, k, crop = NULL){
  if(method == 1){ 
    return(NULL)
  } else {
    tot <- d[reg != k][["tot"]]
    if(method == 2){
      return(tot/max(tot)^0.5)
    } else {
      dqi <- d[reg != k][[paste0(crop, "_dqi")]]  
      return(tot/max(tot)^0.5 * dqi)  
    } 
  }
}
