#' @title getParamLimits
#' @author Oyvind Bleka
#' @description Provides model parameter limitations
#' @details Used further for numerical integration
#' @param mlefit Fitted object using calcMLE
#' @param dev Scaling relative to 2*SE (i.e. 97.5percentil). Should cover most of posterior distr outcome
#' @return returned object from calcINT
#' @export 

getParamLimits = function(mlefit, dev=3) {
  th0 <- mlefit$fit$thetahat #estimate params
  Sigma0 <- mlefit$fit$thetaSigma #estimated covariance of param estimtates
  SE = sqrt(diag(Sigma0)) #obtain standard errors
  lower = th0 - dev*2*SE
  upper = th0 + dev*2*SE
  lower[lower<0] = 0 #Restrict
  scale = abs(mlefit$fit$loglik)
  return(list(lower=lower,upper=upper,scale=scale,acc = NULL))
} #end function

