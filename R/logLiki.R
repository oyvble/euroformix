#' @title logLiki
#' @author Oyvind Bleka
#' @description logLiki returns the likelihood for each marker based on the MLE fit
#'
#' @param mlefit Fitted object using contLikMLE
#' @param verbose Whether printing deconvolution progress. Default is TRUE.
#' @param maxThreads Maximum number of threads to be executed by the parallelization
#' @return ret A vector with log-likelihood-values for each locus for given model
#' @export

logLiki <- function(mlefit,verbose=FALSE,maxThreads=0){  
  return(mlefit$fit$logliki) #already calculated...
}



