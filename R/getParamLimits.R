#' @title getParamLimits
#' @author Oyvind Bleka
#' @description Provides model parameter limitations
#' @details Used further for numerical integration
#' @param mlefit Fitted object using calcMLE
#' @param dev Scaling relative to 2*SE (i.e. 97.5percentil). Should cover most of posterior distr outcome
#' @return returned object from calcINT
#' @export 

getParamLimits = function(mlefit, dev=2) {
  th0 <- mlefit$fit$thetahat #estimate params
  Sigma0 <- mlefit$fit$thetaSigma #estimated covariance of param estimtates
  SE = sqrt(diag(Sigma0)) #obtain standard errors
  if(any(is.na(SE))) stop("Cannot obtain parameter limits since the covariance matrix is not valid. Please optimize a model with less complexity.")
  lower = th0 - dev*2*SE
  upper = th0 + dev*2*SE
  lower[lower<0] = 0 #Restrict
  scale = mlefit$fit$loglik #dont take absolute value
  return(list(lower=lower,upper=upper,scale=scale))
} #end function

