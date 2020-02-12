
#' @title contLikMCMC
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description contLikMCMC simulates from the posterior distribution for a bayesian STR DNA mixture model.
#' @details The procedure are doing MCMC to approximate the marginal probability over noisance parameters. Mixture proportions have flat prior.
#' 
#' If no initial values or covariance matrix has been provided to the function, a call to the MLE function is applied.
#' The Metropolis Hastings routine uses a Multivariate Normal distribution with mean 0 and covariance as delta multiplied with the inverse negative hessian with MLE inserted as transistion kernel.
#' Function calls procedure in c++ by using the package Armadillo and Boost.
#' Marginalized likelihood is estimated using Metropolis Hastings with the "Gelfand and Dey" method.
#'
#' @param mlefit Fitted object using contLikMLE
#' @param uppermu Upper boundary of the mu-parameters
#' @param uppersigma Upper boundary of the sigma-parameters
#' @param upperxi Upper boundary of the xi-parameters
#' @param niter Number of samples in the MCMC-sampling.
#' @param delta A numerical parameter to scale with the covariance function Sigma. Default is 2. Should be higher to obtain lower acception rate.
#' @return ret A list (margL,posttheta,postlogL,logpX,accrat,Ubound ) where margL is Marginalized likelihood for hypothesis (model) given observed evidence, posttheta is the posterior samples from a MC routine, postlogL is sampled log-likelihood values, accrat is ratio of accepted samples. Ubound is upper boundary of parameters.
#' @export 
#' @references Craiu,R.V. and Rosenthal, J.S. (2014). Bayesian Computation Via Markov Chain Monte Carlo. Annu. Rev. Stat. Appl., 1,179-201.
#' @keywords continuous, BayesianModels, MCMC, MetropolisHastings, MarginalizedLikelihoodEstimation

contLikMCMC = function(mlefit,niter=1e4,delta=10,maxxi=1) {
 #A mlefit object returned from contLikMLE is required to do MCMC!
 loglik0 <-  mlefit$fit$loglik #get maximized likelihood
 model <- mlefit$model
 th0 <- mlefit$fit$thetahat
 Sigma0 <- mlefit$fit$thetaSigma
 varnames <- names(mlefit$fit$thetahat) #variable names
 if(!all(length(th0)%in%dim(Sigma0))) stop("Length of th0 and dimension of Sigma was not the same!")
 ret <- prepareC(model$nC,model$samples,model$popFreq,model$refData,model$condOrder,model$knownRef,model$kit)
 nC <- ret$nC
 np <- length(th0) #number of unknown parameters
 nodeg <- is.null(model$kit) #check for degradation

 if(is.null(model$xi)) {
   loglikYtheta <- function(theta) {   #call c++- function: length(phi)=nC+1
    if(any(theta<0)) return(-Inf) 
    xi1 <- theta[np] #value of xi
    if(theta[np]>maxxi) return(-Inf) #special case for xi which has a upper boundary
    if(nodeg) theta <- c(theta[1:(nC+1)],1,xi1)
    Cval  <- .C("loglikgammaC",as.numeric(0),as.numeric(theta),as.integer(np),ret$nC,ret$nK,ret$nL,ret$nS,ret$nA,ret$obsY,ret$obsA,ret$CnA,ret$allAbpind,ret$nAall,ret$CnAall,ret$Gvec,ret$nG,ret$CnG,ret$CnG2,ret$pG,ret$pA, as.numeric(model$prC), ret$condRef,as.numeric(model$threshT),as.numeric(model$fst),ret$mkvec,ret$nkval,as.numeric(model$lambda),ret$bp,as.integer(0),PACKAGE="euroformix")[[1]]
    loglik <- Cval + model$pXi(xi1) #weight with prior of tau and 
    return(loglik) #weight with prior of tau and stutter.
   }
 } else {  
   loglikYtheta <- function(theta2) {   #call c++- function: length(phi)=nC
    if(nodeg) theta2 <- c(theta2,1)
    theta <- c(theta2,model$xi) #stutter-parameter added as known
    if(any(theta<0)) return(-Inf) 
    Cval  <- .C("loglikgammaC",as.numeric(0),as.numeric(theta),as.integer(np),ret$nC,ret$nK,ret$nL,ret$nS,ret$nA,ret$obsY,ret$obsA,ret$CnA,ret$allAbpind,ret$nAall,ret$CnAall,ret$Gvec,ret$nG,ret$CnG,ret$CnG2,ret$pG,ret$pA, as.numeric(model$prC), ret$condRef,as.numeric(model$threshT),as.numeric(model$fst),ret$mkvec,ret$nkval,as.numeric(model$lambda),ret$bp,as.integer(0),PACKAGE="euroformix")[[1]]
    return(Cval)
   }
 }
 C <- chol(delta*Sigma0) #scale variance with a factor 2: ensures broad posterior
 X <- t( t(C)%*%matrix(rnorm(np*niter),ncol=niter,nrow=np)) #proposal values

 logdmvnorm <- function(X,mean,cholC) { #function taken from mvtnorm-package
   p <- nrow(cholC)
   tmp <- backsolve(cholC,t(X)-mean,p,transpose=TRUE)
   rss <- colSums(tmp^2)
   logretval <- -sum(log(diag(cholC))) - 0.5*p*log(2*pi) - 0.5*rss
   return(logretval)
 }
 #removed:  Importance sampling using Normal(th0,delta*Sigma)
 if(1) { #MCMC by Gelfand and Dey (1994), using h() = Normal(th0,delta*Sigma)
   rlist <- list()
   if(nC>1) rlist[[length(rlist)+1]] <- 1:(nC-1)
   rlist[[length(rlist)+1]] <- nC:(nC+1+!nodeg)
   if(is.null(model$xi))  rlist[[length(rlist)+1]] <- np
   nB <- length(rlist) #number of blocks
   M2 <- nB*niter+1
   postth <- matrix(NA,ncol=np,nrow=M2) #accepted th
   postlogL <- rep(NA,M2) #accepted th
   postth[1,] <- th0
   postlogL[1] <- loglik0 #loglikYtheta(th0) #get start-likelihood   
   U <- runif(M2) #random numbers
   m <- 2 #counter for samples
   m2 <- 1 #counter for proposal
   nacc <- 0  
   while(m<=M2) {
    for(r in 1:nB ) { #for each blocks
     range <- rlist[[r]]
     postth[m,] <-  postth[m-1,] #proposed th
     postth[m,range] <- X[m2,range] + postth[m,range]
     postlogL[m] <- loglikYtheta(postth[m,])
     pr <- exp(postlogL[m]- postlogL[m-1]) #acceptance rate
     if(U[m]>pr) { #if not accepted, i.e. random pr too large
      postth[m,] <-  postth[m-1,]
      postlogL[m] <- postlogL[m-1 ]
     } else {
      nacc <- nacc + 1
     }
     m <- m + 1 #update counter
    } #end for each blocks
    m2 <- m2 +1 #update proposal counter
   } #end while not done
  accrat <- nacc/M2 #acceptance ratio
  logpX <- logdmvnorm(postth,mean=th0,cholC=chol(Sigma0)) #insert with Normal-approx of post-th
  #plot(postlogL,ty="l")
  #plot(logpX,ty="l")
  margL <- 1/mean(exp(logpX - postlogL)) #estimated marginal likelihood
 }
# nU <- nC-ret$nK #number of unknowns
# if(nU>1) { #if more than 1 unknown 
#  margL <- factorial(nU)*margL #get correct ML adjusting for symmetry
# } #end method
 colnames(postth) <- varnames  #save variable names
  return(list(margL=margL,posttheta=postth,postlogL=postlogL,logpX=logpX,accrat=accrat,MLE=mlefit$fit$thetahat,Sigma=mlefit$fit$thetaSigma))
} #end function

