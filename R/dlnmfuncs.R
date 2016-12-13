#' Stratified lambda
#' 
#' This is a custom lambda function to enter into eesim.  It creates lambda 
#' values based on stratified relative risks.  
#' 
#' @examples 
#' test <- rep(0, 30)
#' test[c(12,14,23)]<-1
#' testbasis <- crossbasis(test, lag = 2, argvar = list(fun="lin"), 
#'            arglag = list(fun="strata", breaks = c(0,1,2)))[,]
#' logrr <- testbasis[,2:ncol(testbasis)]%*%rrmat
#' stratalambda(test, baseline = rep(58, 30), rrvalues = c(1.5,1.2,1.05), 
#'              lag = 2, argvar = list(fun="lin"), 
#'              arglag = list(fun="strata", breaks = c(0,1,2)))
#' 
#' @export
#' 
stratalambda <- function(exposure, baseline, rrvalues, lag, argvar, arglag){
  #Make basis of ones where there is any relative risk
  basis <- crossbasis(exposure, lag=lag, argvar=argvar,arglag=arglag)
  #Turn basis into matrix
  basis <- as.matrix(basis)
  #Take log of relative risk values, make rx1 matrix from them:
  rrmat <- matrix(log(rrvalues), ncol = 1)
  #Matrix multiplication: results in nx1 vector of log(rr), adds overlapping log(rr):
  logrr <- basis[ , 2:ncol(basis)] %*% rrmat
  #Create vector zeros:
  log_lambda <- rep(0, length(logrr))
  #Something in here is incorrect:
  for (i in 1:length(logrr)){
    if (is.na(logrr[i])==T){
      log_lambda[i] <- log(baseline[i])
    }
    else if (logrr[i] != 0){
      log_lambda[i] <- log(baseline[i]) + logrr[i]
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
#' testexpo <- sim_exposure(n=200, central = .01, amp = .01, 
#'             exposure_type = "binary")
#' testrr <- smoothrr(testexpo$x, lag = 20, scale = 6)
#' overlap <- rep(0,200)
#' overlap[c(50, 60, 144)] <-1
#' smoothrr(testexpo$x, lag = 20)
#' 
#' @export
#' 
smoothrr <- function(exposure, lag, scale=6){
  #scale of 6 results in max rr of 1.18, independent of other arguments.
  #Add lag amount of burn time at beginning of exposure
  exposure <- c(rep(0, lag), exposure)
  #Create a vector that has a 1 on each day there was increased rr:
  lagtime <- dlnm::crossbasis(x=exposure, lag = lag, argvar = list(fun = "thr", 
                              thr.value = 0))[,]
  rrseq <- rep(0, length(lagtime))
  for (i in 1:length(lagtime)){
    if(is.na(lagtime[i])){
      rrseq[i] <- 0
    }
  else if (lagtime[i]==1 & (lagtime[i-1]==0 | is.na(lagtime[i-1]))){
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