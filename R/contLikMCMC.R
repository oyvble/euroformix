
#' @title contLikMCMC
#' @author Oyvind Bleka
#' @description calcMCMC provides samples from the posterior distribution for the model.
#' @details The procedure are doing MCMC to approximate the marginal probability over model parameters. 
#' 
#' The Metropolis Hastings routine uses following proposal: Multivariate Normal distribution with mean 0 and covariance as delta multiplied with the inverse negative hessian with MLE inserted.
#' Marginalized likelihood (Bayesian) is estimated using Metropolis Hastings with the "GD-method, Gelfand and Dey (1994).
#'
#' @param mlefit Fitted object using calcMLE
#' @param niter Number of samples in the MCMC-sampling.
#' @param delta A numerical parameter to scale with the covariance function Sigma. Default is 2. Should be higher to obtain lower acception rate.
#' @param verbose Whether printing simulation progress. Default is TRUE
#' @param seed The user can set seed if wanted
#' @param usePhi Whether the transformed domain of the parameters (phi) should be used instead. This affects the priors!
#' @param mcmcObj An object from contLikMCMC output
#' @return ret A list (logmargL,posttheta,postlogL,logpX,accrat) where logmargL is Marginalized likelihood, posttheta is the posterior samples, postlogL is log-likelihood values
#' @export 

#niter=1e4;delta=2;verbose=TRUE;seed=NULL;usePhi=FALSE;mcmcObj=NULL
contLikMCMC = function(mlefit,niter=1e4,delta=2,verbose=TRUE,seed=NULL,usePhi=FALSE, mcmcObj=NULL) {
  if(!is.null(seed)) {
    if(!is.null(mcmcObj)) seed = seed + length(mcmcObj$postlogL) #adjust seed wrt number of samples
    set.seed(seed) #set seed if provided    
  } 
  loglik0 <-  mlefit$fit$loglik #get maximized likelihood
  model <- mlefit$model
  priorBWS = model$priorBWS #obtain prior stutter model
  priorFWS = model$priorFWS #obtain prior stutter model
  modTypes <- mlefit$modelType #obtain model types
  
  th0 <- mlefit$fit$thetahat #estimate params
  Sigma0 <- mlefit$fit$thetaSigma #estimated covariance of param estimtates
  if(usePhi) {
    th0 <- mlefit$fit$phihat #MLE
    Sigma0 <- mlefit$fit$phiSigma #estimated covariance of param estimtates
  }
  np <- length(th0) #number of unknown parameters
  if(np!=ncol(Sigma0)) stop("Length of th0 and dimension of Sigma was not the same!")
  c <- mlefit$prepareC #returned from prepareC
  nC = model$nC #number of contributors
  
  #PREPARE:
  if(verbose) print("Carrying out preparations for MCMC")
  mod = Rcpp::Module( "mod",PACKAGE="euroformix" ) #load module
  obj = methods::new(mod$ExposedClass) #create object of class
  
  #Step 1: insert data to exposed class (filldata)
  obj$filldata(c$nStutterModels,c$nMarkers,c$nRepMarkers,c$nAlleles,c$startIndMarker_nAlleles,c$startIndMarker_nAllelesReps,c$peaks,c$freqs,c$dropinWeight, c$nTyped, c$maTyped, c$basepair,
              c$BWfrom, c$FWfrom, c$BWto, c$FWto, c$nPotStutters, c$startIndMarker_nAllelesTot, c$QalleleIndex, c$dropinProb, c$fst, c$AT, c$NOK, c$knownGind, c$relGind, c$ibd, mlefit$maxThreads) 
  
  #Step 2: Indexing large matrix (doIndex)
  obj$prepare(as.integer(nC))
  
  #The likelihood function taking theta-parameters (non-transformed params)
  loglikphi <- function(par) {
    #   phi=postth[m,]
    param = .convBack(par,nC, modTypes,isPhi=usePhi) #need to convert params back
    if(any(param<0) || any(param[-(nC+1:2)]>1)) return(-Inf) #outside of domain   
    loglik = obj$loglik(as.numeric(param) )
    if(!is.null(priorBWS) && param[nC+4]>0)  loglik <- loglik + log(priorBWS(param[nC+4])) #weight with prior of xi
    if(!is.null(priorFWS) && param[nC+5]>0)  loglik <- loglik + log(priorFWS(param[nC+5])) #weight with prior of xiFW
    if(is.nan(loglik)) return(-Inf) #avoid NAN values
    
    if(verbose) { #only show progressbar if verbose
     progcount <<- progcount + 1
     setTxtProgressBar(progbar,progcount)
    } 
    return(loglik) #weight with prior of tau and stutter.
  }
  
  if(!is.null(mcmcObj)) delta = mcmcObj$delta #important that same delta is used
  C <- chol(delta*Sigma0) #scale variance with a factor 2: ensures broad posterior
  X <- matrix(rnorm(np*niter),ncol=np,nrow=niter)%*%C #proposal values
  
  logdmvnorm <- function(X2,mean,cholC) { #function taken from mvtnorm-package
    p <- nrow(cholC)
    tmp <- backsolve(cholC,X2-mean,p,transpose=TRUE)
    rss <- colSums(tmp^2)
    logretval <- -sum(log(diag(cholC))) - 0.5*p*log(2*pi) - 0.5*rss
    return(logretval)
  }
  
  #Inititate progressbar if verbose
  progcount = 1  #counter
  if( verbose ) progbar <- txtProgressBar(min = 0, max = niter, style = 3) #create progress bar
  
  #Performing MCMC Metropolis Hastings by Gelfand and Dey (1994), using h() = Normal(th0,delta*Sigma)
  postth <- matrix(NA,ncol=np,nrow=niter) #accepted th
  postlogL <- rep(NA,niter) #accepted th
  if(!is.null(mcmcObj)) { #if we continue the same chain (use info from previous run)
    niterPrev = length(mcmcObj$postlogL) #number of prev iterations
    postlogL[1] = mcmcObj$postlogL[niterPrev]
    th1 = mcmcObj$posttheta[niterPrev,-nC] #obtain previously samples
    postth[1,] = th1[1:np] # obtain last accepted sample (last state)
  } else { #if this is first run
    th1 = th0 #start adjustment
    postth[1,] <- X[1,] + th0 #add noise into start point: IMPORTANT TO START VARIATION EARLIER
    postlogL[1] <- loglikphi(postth[1,])# loglik0 #init with proper value logliktheta(postth[1,]) #get start-likelihood  
  }
  if(is.infinite(postlogL[1])) { #fall back to optimum if failed (unlikely if model is proper)
    postth[1,] = th0
    postlogL[1] = loglik0
  }
  U <- runif(niter) #random numbers
  m <- 2 #counter for samples
  nacc <- 0  
  while(m<=niter) {
    postth[m,] <-  X[m,] + postth[m-1,] #proposed theta
    postlogL[m] <- loglikphi(postth[m,])
    pr <- exp(postlogL[m]- postlogL[m-1]) #acceptance rate
    if(U[m]>pr) { #if not accepted, i.e. random prob too large (above pr)
      postth[m,] <-  postth[m-1,]
      postlogL[m] <- postlogL[m-1 ]
    } else {
      nacc <- nacc + 1
    }
    m <- m + 1 #update counter
  } #end while not done
  obj$close() #free memory
  if( verbose ) cat("\n") #new line
  
  #POSTE PROCESSING:
  accrat <- nacc/(niter-1) #acceptance ratio (dont count first)
  logpX <- logdmvnorm(t(postth),mean=th0,cholC=C) #insert with Normal-approx of post-th
  if(!is.null(mcmcObj)) { #ADDING RESULTS FROM EARLIER RUN
    postlogL = c(mcmcObj$postlogL,postlogL) #append vector
    logpX = c(mcmcObj$logpX,logpX) #append vector
  }
  
  logVals = logpX - postlogL #this is logged values intended to be "exped"
  #plot(postlogL,ty="l")
  #plot(logpX,ty="l")
  offset = max(logVals) #find max value
  #margL <- 1/mean(exp(logVals - offset)) #estimated marginal likelihood
  #logmargL = log(margL) - offset #this is identical as the next line:
  logmargL <- log(length(logVals)) - offset - log(sum(exp(logVals-offset))) #estimated marginal likelihood (logged)
  
  #convert back proposed values to theta-format
  postth = .convBack(postth, nC,modTypes,isPhi = usePhi)
  if(!is.null(mcmcObj))  postth = rbind(mcmcObj$posttheta,postth) #append matrix

  #Whether to calculate sub-source LR (permutating the position of contributor):  
  # nU <- nC-ret$nK #number of unknowns
  # if(nU>1) { #if more than 1 unknown 
  #  margL <- factorial(nU)*margL #get correct ML adjusting for symmetry
  # } #end method
  colnames(postth) <-  names(mlefit$fit$thetahat2)  #save variable names
  return(list(logmargL=logmargL,posttheta=postth,postlogL=postlogL,logpX=logpX,accrat=accrat,MLE=mlefit$fit$thetahat2,Sigma=mlefit$fit$thetaSigma2,
             offset=offset,seed=seed,delta=delta,usePhi=usePhi))
} #end function

