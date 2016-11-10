#' Stratified lambda
#' 
#' This is a custom lambda function to enter into eesim.  It creates lambda 
#' values based on stratified relative risks.  
#' 
#' 
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
#' Continuous Relative Risk
#' 
#' This is an intermediate function to use while making the smooth rr custom 
#' lambda function for eesim.
#' 
#' @examples 
#' set.seed(15)
#' testexpo <- sim_exposure(n=200, central = .01, amp = .01, exposure_type = "binary")
#' testrr <- smoothrr(testexpo$x, lag = 20, scale = 6)
#' 
#' @export
#' 
smoothrr <- function(exposure, lag, scale){
  #Create a vector that has a 1 on each day there was increased rr:
  lagtime <- dlnm::crossbasis(x=exposure, lag = lag, argvar = list(fun = "thr", thr.value = 0))[,]
  rrseq <- rep(0, length(lagtime))
  for (i in 1:length(lagtime)){
    if(is.na(lagtime[i])){
      rrseq[i] <- 0
    }
  else if (lagtime[i]==1 & lagtime[i-1]==0){
    rrseq[i:(i+lag-1)] = abs(-lag:-1)
  }
  else if (lagtime[i]==0){
    rrseq[i]<-0
  }
}
polylag <- dlnm:::poly(rrseq)/scale
rr <- exp(polylag)
return(rr)
}
#' 
#' Smooth lambda
#' 
#' This function can be input as a custom lambda function in eesim.  It uses the
#' smoothrr function to create lambda values based on smooth relative risks.
#' 
#' @examples 
#' set.seed(15)
#' testexpo <- sim_exposure(n=200, central = .01, amp = .01, exposure_type = "binary")
#' base <- create_baseline(200, 58, trend = "no trend")
#' smoothlambda(exposure=testexpo$x, baseline = base$baseline, lag = 20, scale = 2)
#' 
#' @export
#' 
smoothlambda <- function(exposure, baseline, lag, scale){
  rr <- smoothrr(exposure, lag, scale)
  log_lambda <- log(baseline)+log(rr)
  lambda <- exp(log_lambda)
  return(lambda)
}
#' 
#' 
#' 
#' 
#' 