# measures of accuracy functions
# x or x[,1] is observed, y or x[,2] is predicted

.xy2dt <- function(x,y){
  err.msg <- "input should be a 2 column matrix/df/DT or 2 vectors"
  if(is.null(y)){
    if(!identical(ncol(x), 2L)) stop(err.msg)
    if(is.matrix(x)) x <- as.data.frame(x)
    setDT(x)
    setnames(x, names(x), c("x", "y"))
    return(x)
  } else {
    if(!all(sapply(list(x,y), is.vector, mode = "numeric"))) stop(err.msg)
    return(data.table(x = unname(x), y = unname(y)))
  }
}

f_mse <- function(x, y=NULL, na.rm=TRUE){
  dt <- .xy2dt(x=x, y=y)
  dt[, mean((x - y)^2, na.rm = na.rm)]
}

f_mse_ov <- function(x, y=NULL, obs_1st=TRUE, na.rm=TRUE){
  dt <- .xy2dt(x=x, y=y)
  if(obs_1st){
    dt[, mean((pmin(x - y, 0))^2, na.rm = na.rm)]  
  } else {
    dt[, mean((pmin(y - x, 0))^2, na.rm = na.rm)]  
  }
}

f_mse_un <- function(x, y=NULL, obs_1st=TRUE, na.rm=TRUE){
  dt <- .xy2dt(x=x, y=y)
  if(obs_1st){
    dt[, mean((pmax(x - y, 0))^2, na.rm = na.rm)]  
  } else {
    dt[, mean((pmax(y - x, 0))^2, na.rm = na.rm)]  
  }
}

f_mse_pr <- function(x, y=NULL, obs_1st=TRUE, pres_th=1e-5, na.rm=TRUE){
  dt <- .xy2dt(x=x, y=y)
  dt[, `:=`(xpr = x, ypr = y)]
  if(obs_1st){
    dt[x < pres_th, `:=`(xpr = 0, ypr = 0)]
  } else {
    dt[y < pres_th, `:=`(xpr = 0, ypr = 0)]
  }
  out <- dt[, mean((xpr - ypr)^2, na.rm = na.rm)]
  dt[,c("xpr", "ypr"):= NULL]
  return(out)
}

f_mse_ab <- function(x, y=NULL, obs_1st=TRUE, pres_th=1e-5, na.rm=TRUE){
  dt <- .xy2dt(x=x, y=y)
  dt[,`:=`(xab = x, yab = y)]
  if(obs_1st){
    dt[x >= pres_th, `:=`(xab = 0, yab = 0)]
  } else {
    dt[y >= pres_th, `:=`(xab = 0, yab = 0)]
  }
  out <- dt[, mean((xab - yab)^2, na.rm = na.rm)]
  dt[,c("xab", "yab"):= NULL]
  return(out)
}

f_mofa <- function(x, y=NULL, obs_1st=TRUE, pres_th=1e-5, na.rm=TRUE){
  dt <- .xy2dt(x=x, y=y)
  mse <- f_mse(dt)
  mse_ov <- f_mse_ov(dt, obs_1st=obs_1st, na.rm=na.rm)
  mse_un <- f_mse_un(dt, obs_1st=obs_1st, na.rm=na.rm)
  mse_ab <- f_mse_ab(dt, obs_1st=obs_1st, pres_th=pres_th, na.rm=na.rm)
  mse_pr <- f_mse_pr(dt, obs_1st=obs_1st, pres_th=pres_th, na.rm=na.rm)
  data.table(stat = c("mse", "mse_ov", "mse_un", "mse_ab", "mse_pr"),
             value = c(mse, mse_ov, mse_un, mse_ab, mse_pr))
             
}

