#' Stratified lambda
#' 
#' This is a custom lambda function to enter into eesim.  It creates lambda 
#' values based on stratified relative risks.  
#' 
#' @export
#' 
stratalambda <- function(exposure, baseline, rrvalues, lag, argvar, arglag){
  basis <- crossbasis(exposure, lag=lag, argvar=argvar,arglag=arglag)
  basis <- as.data.frame(basis)
  for (i in 2:ncol(basis)){
    basis[,i] <- basis[,i]*rrvalues[i-1]
  }
  rr <- apply(basis, MARGIN = 1, FUN = sum, na.rm = F)
  log_lambda <- rep(0, length(rr))
  for (i in 1:length(rr)){
    if (is.na(rr[i])==T){
      log_lambda[i] <- log(baseline[i])
    }
    else if (rr[i] != 0){
      log_lambda[i] <- log(baseline[i]) + log(rr[i])
    }
    else{
      log_lambda[i] <- log(baseline[i])
    }
  }
  lambda <- exp(log_lambda)
  return(lambda)
}

#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 
#' 