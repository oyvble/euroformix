
#' @title contLikMCMCpara
#' @author Oyvind Bleka
#' @description Same as the contLikMCMC function (dummy function)
#' 
#' @param mlefit Fitted object using contLikMLE
#' @param niter Number of samples in the MCMC-sampling.
#' @param delta A numerical parameter to scale with the covariance function Sigma. Default is 2. Should be higher to obtain lower acception rate.
#' @param maxxi Upper boundary of the xi-parameters
#' @param maxxiFW Upper boundary of the xiFW-parameters
#' @param verbose Boolean whether printing simulation progress. Default is TRUE
#' @param seed The user can set seed if wanted
#' @param maxThreads Maximum number of threads to be executed by the parallelization
#' @return ret A list (margL,posttheta,postlogL,logpX,accrat,Ubound ) where margL is Marginalized likelihood for hypothesis (model) given observed evidence, posttheta is the posterior samples from a MC routine, postlogL is sampled log-likelihood values, accrat is ratio of accepted samples. Ubound is upper boundary of parameters.
#' @export 

contLikMCMCpara = function(mlefit,niter=1e4,delta=2,maxxi=1,maxxiFW=1,verbose=TRUE,seed=NULL,maxThreads=32) {

  return( contLikMCMC(mlefit=mlefit,niter=niter,delta=delta,maxxi=maxxi,maxxiFW=maxxiFW,verbose=verbose,seed=seed,maxThreads=maxThreads) )
  
} #end function

