#' @title contLikMCMC
#' @author Oyvind Bleka
#' @description calcMCMC provides samples from the posterior distribution of the model parameters
#' @details The procedure also uses the samples to approximate the marginal likelihood used for Bayes Factor 
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
  thetaSE = sqrt(diag(Sigma0)) #obtain diagnoals (variance)
  #indicate which params should be fixed:  
  isFixed = thetaSE < 1e-3 #Criterion similar as using dev=~2 in integrals (0.01 width)
  
  #Obtain cholesky decomposition of covariance matrix
  thetaFixed = th0[isFixed] #fixate these params
  thetaVary = th0[!isFixed] #theta params to vary
  Sigma0 = Sigma0[!isFixed,!isFixed,drop=FALSE]
  C = NULL
  tryCatch({
    C <- chol(delta*Sigma0) #scale variance with a factor 2: ensures broad posterior
  },error=function(e) e)
  if(is.null(C)) {
    if(verbose) print("NOTE: Covariance matrix was modified to enable MCMC runs...") 
    #diag(Sigma0)
    #if(is.null(C))  stop("MCMC could not be carried out since cholesky decomposition of the covariance matrix failed. Please fit a model with reduced complexity and try again!")
    Sigma0 <- Sigma0 + 0.01 #small values added to covariance matrix
    C <- chol(delta*Sigma0) #scale variance with a factor 2: ensures broad posterior
  } 

  np <- length(thetaVary) #number of unknown parameters (to vary)
  if(np!=ncol(Sigma0)) stop("Length of th0 and dimension of Sigma was not the same!")
  c <- mlefit$prepareC #returned from prepareC
  nC = model$nC #number of contributors

  #PREPARE:
  if(verbose) print("Carrying out preparations for MCMC")
  obj = prepareCobj(c, verbose) #use wrapper function to obtain C++ pointer
  
  restTime=system.time({ #prepare restriction (can take a while)
    preSearched = .preSearch(obj,c,mlefit$thetaPresearch,mlefit$resttol,priorBWS, priorFWS)
  })[3]
  if(verbose) print(paste0("Presearch done and took ",round(restTime),"s. Start simulate..."))

  #The likelihood function taking theta-parameters (non-transformed params)
  logLik <- function(parVary) {
    #   phi=postth[m,]
    par = rep(NA,length(th0)) #init param vector 
    par[!isFixed] = parVary #insert varying values
    par[isFixed] = thetaFixed #insert fixed values
    param = .convBack(par,nC, modTypes,isPhi=usePhi) #Obtain full vector of params
    if(any(param<0) || any(param[-(nC+1:2)]>1)) return(-Inf) #outside of domain   
    mxValid = .checkMxValid(param[1:nC],c$nC,c$nU,c$hasKinship)
    if(!mxValid) return(-Inf)

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
    postth[1,] <- X[1,] + thetaVary #add noise into start point: IMPORTANT TO START VARIATION EARLIER
    postlogL[1] <- logLik(postth[1,])# loglik0 #init with proper value logliktheta(postth[1,]) #get start-likelihood  
  }
  if(is.infinite(postlogL[1])) { #fall back to optimum if failed (unlikely if model is proper)
    postth[1,] = thetaVary
    postlogL[1] = loglik0
  }
  U <- runif(niter) #random numbers
  m <- 2 #counter for samples
  nacc <- 0  
  while(m<=niter) {
    postth[m,] <-  X[m,] + postth[m-1,] #proposed theta
    postlogL[m] <- logLik(postth[m,])
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
  logpX <- logdmvnorm(t(postth),mean=thetaVary,cholC=C) #insert with Normal-approx of post-th
  if(!is.null(mcmcObj)) { #ADDING RESULTS FROM EARLIER RUN
    postlogL = c(mcmcObj$postlogL,postlogL) #append vector
    logpX = c(mcmcObj$logpX,logpX) #append vector
  }
  
  logVals = logpX - postlogL #this is logged values intended to be "exped"
# plot(postlogL,ty="l")
# plot(logpX,ty="l")
  offset = max(logVals) #find max value and use this as offset
  #margL <- 1/mean(exp(logVals - offset)) #estimated marginal likelihood
  #logmargL = log(margL) - offset #this is identical as the next line:
  logmargL <- log(length(logVals)) - offset - log(sum(exp(logVals-offset))) #estimated marginal likelihood (logged)
  logCombAdjFactor = 0 #no adjustment
  nUnknown = c$nU - as.integer(c$hasKinship) #number of Mx params to be restricted
  if(nUnknown>1) logCombAdjFactor = lfactorial(nUnknown)   #Need to take into account the symmetry of unknowns  
  logmargL = logmargL + logCombAdjFactor #make adjustment
  #convert back proposed values to theta-format
  #Special handling if there was a variable to fix: need to put back
  if(any(isFixed)) {
    postthAll = matrix(NA,ncol=niter,nrow=length(th0)) #must be this direction to ease the fill in
    postthAll[isFixed,] = thetaFixed  #put in f
    postthAll = t(postthAll) #transpose back
    postthAll[,!isFixed] = postth #insert other params
    postth = .convBack(postthAll, nC,modTypes,isPhi = usePhi)
  } else {
    postth = .convBack(postth, nC,modTypes,isPhi = usePhi)
  }

  if(!is.null(mcmcObj))  postth = rbind(mcmcObj$posttheta,postth) #append matrix

  colnames(postth) <-  names(mlefit$fit$thetahat2)  #save variable names
  return(list(logmargL=logmargL,posttheta=postth,postlogL=postlogL,logpX=logpX,accrat=accrat,MLE=mlefit$fit$thetahat2,Sigma=mlefit$fit$thetaSigma2,
             offset=offset,seed=seed,delta=delta,usePhi=usePhi,logCombAdjFactor=logCombAdjFactor))
} #end function

