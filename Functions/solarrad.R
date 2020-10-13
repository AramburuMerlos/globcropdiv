# function derived from https://github.com/cran/sirad/tree/master/R
rta <- function(lat_degrees, month){
  if(is.na(lat_degrees)) return(NA)
  
  Con <- 4.921
  days_per_month <- c(31,28.25,31,30,31,30,31,31,30,31,30,31) 
  doy_in_month <- cumsum(round(days_per_month))
  doys <- ifelse(month == 1,1,(doy_in_month[month - 1] + 1)):doy_in_month[month]
  
  lat_rad <- lat_degrees * pi/180
  
  # Corrected Earth Sun Distance
  cesd <- 1 + 0.0334 * cos(0.01721 * doys - 0.0552) 
  # Solar Decline
  soldec <- 0.4093 * sin((2 * pi * (284 + doys))/365)
  if (abs(lat_degrees) < 66.5) {
    # Daylight Time Factor
    dltm <- acos(-tan(lat_rad) * tan(soldec))
    Sd <- Con*24/pi*cesd*(sin(lat_rad)*sin(soldec)*dltm + cos(lat_rad)*cos(soldec)*sin(dltm))
  } else { 
    # Solar zenith angle (hourly)
    sza <- matrix(ncol = length(doys), nrow = 24)
    tr <- 0.2618
    for(hr in 1:24){
      sza[hr,] <- acos(sin(lat_rad)*sin(soldec) + cos(lat_rad)*cos(soldec)*cos(tr*(hr-12)))
    }
    shr <- apply(sza, 1, function(x) Con * cesd * cos(x))
    shr[shr<0] <- 0      
    Sd <- apply(shr,1, sum) 
  }
  Sd[Sd<0] <- 0
  mean(Sd)  
}



