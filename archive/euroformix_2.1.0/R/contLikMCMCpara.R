
#' @title contLikMCMCpara
#' @author Oyvind Bleka
#' @description Parallelization on contLikMCMC using snow
#' @details The procedure is doing parallelization of the contLikMCMC function
#' 
#' @param mlefit Fitted object using contLikMLE
#' @param niter Number of samples in the MCMC-sampling.
#' @param delta A numerical parameter to scale with the covariance function Sigma. Default is 2. Should be higher to obtain lower acception rate.
#' @param maxxi Upper boundary of the xi-parameters
#' @return ret A list (margL,posttheta,postlogL,logpX,accrat,Ubound ) where margL is Marginalized likelihood for hypothesis (model) given observed evidence, posttheta is the posterior samples from a MC routine, postlogL is sampled log-likelihood values, accrat is ratio of accepted samples. Ubound is upper boundary of parameters.
#' @export 
#' @references Craiu,R.V. and Rosenthal, J.S. (2014). Bayesian Computation Via Markov Chain Monte Carlo. Annu. Rev. Stat. Appl., 1,179-201.
#' @keywords continuous, BayesianModels, MCMC, MetropolisHastings, MarginalizedLikelihoodEstimation

contLikMCMCpara = function(mlefit,niter=1e4,delta=10,maxxi=1) {
 library(parallel, warn.conflicts = FALSE)

 ncores <- parallel::detectCores() #number of physical cores (parallel)=number of chains
 nCl <- ncores #number of clusters
 print(paste0("\nNumber of paralell chains will be ",nCl ))
 #if( nCl==1 || nCl%%2!=0) stop("Please change the number of startpoints to an even number. This is can be changed under Optimization.")
 cl <- parallel::makeCluster(nCl,type="PSOCK")
 inputlist <- list(mlefit=mlefit,niter=niter,delta=delta,maxxi=maxxi)
 inputlist <- rep(list(inputlist),nCl) #number of clusters

 fclust<-function(x) {
  library(euroformix)
  return(euroformix::contLikMCMC(x$mlefit,x$niter,x$delta,x$maxxi))
 }
 retlist <- parallel::clusterApply(cl,inputlist,fclust)
 parallel::stopCluster(cl)

 #Merge samples. Evt. averaging over end-results
#names(retlist[[1]])
 retlist2 <- retlist[[1]]
 if(nCl>1) { #if at least 1 parallel
  for(c in 2:nCl) {
    retlist2$posttheta <- rbind(retlist2$posttheta,retlist[[c]]$posttheta) 
    retlist2$postlogL <- c(retlist2$postlogL,retlist[[c]]$postlogL) 
    retlist2$logpX <- c(retlist2$logpX,retlist[[c]]$logpX) 
  }
  retlist2$margL <- mean(sapply(retlist,function(x) x$margL))
  retlist2$accrat <- mean(sapply(retlist,function(x) x$accrat))
 }
 return(retlist2)
} #end function

