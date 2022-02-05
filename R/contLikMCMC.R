
#' @title contLikMCMC
#' @author Oyvind Bleka
#' @description contLikMCMC provides samples from the posterior distribution for the model.
#' @details The procedure are doing MCMC to approximate the marginal probability over model parameters. 
#' 
#' The Metropolis Hastings routine uses following proposal: Multivariate Normal distribution with mean 0 and covariance as delta multiplied with the inverse negative hessian with MLE inserted.
#' Marginalized likelihood (Bayesian) is estimated using Metropolis Hastings with the "GD-method, Gelfand and Dey (1994).
#'
#' @param mlefit Fitted object using contLikMLE
#' @param niter Number of samples in the MCMC-sampling.
#' @param delta A numerical parameter to scale with the covariance function Sigma. Default is 2. Should be higher to obtain lower acception rate.
#' @param maxxi Upper boundary of the xi-parameters
#' @param maxxiFW Upper boundary of the xiFW-parameters
#' @param verbose Whether printing simulation progress. Default is TRUE
#' @param seed The user can set seed if wanted
#' @param maxThreads Maximum number of threads to be executed by the parallelization
#' @return ret A list (logmargL,posttheta,postlogL,logpX,accrat,Ubound ) where margL is Marginalized likelihood for hypothesis (model) given observed evidence, posttheta is the posterior samples from a MC routine, postlogL is sampled log-likelihood values, accrat is ratio of accepted samples. Ubound is upper boundary of parameters.
#' @export 
#' @examples
#' \dontrun{
#' kit = "ESX17"
#' sep0 = .Platform$file.sep
#' AT0 = 50
#' popfn = paste(path.package("euroformix"),"FreqDatabases",paste0(kit,"_Norway.csv"),sep=sep0)
#' evidfn = paste(path.package("euroformix"),"examples",paste0(kit,"_3p.csv"),sep=sep0)
#' reffn = paste(path.package("euroformix"),"examples",paste0(kit,"_refs.csv"),sep=sep0)
#' popFreq = freqImport(popfn)[[1]] #obtain list with population frequencies
#' samples = sample_tableToList(tableReader(evidfn))
#' refData = sample_tableToList(tableReader(reffn))
#' dat = prepareData(samples,refData=refData,popFreq=popFreq,threshT=AT) 
#' plotEPG2(dat$samples,dat$refData,kit=kit,AT=AT0)
#' mlefit = contLikMLE(3,dat$samples,dat$popFreq,dat$refData,1:3,
#' 	kit=kit,xi=NULL,prC=0.05,lambda=0.01,seed=1)
#' mcmcfit = contLikMCMC(mlefit,niter=5000,delta=1,seed=1) #mcmcfit$acc
#' }


contLikMCMC = function(mlefit,niter=1e4,delta=2,maxxi=1,maxxiFW=1,verbose=TRUE,seed=NULL,maxThreads=32) {
 if(!is.null(seed)) set.seed(seed) #set seed if provided
  
 loglik0 <-  mlefit$fit$loglik #get maximized likelihood
 model <- mlefit$model
 xi = model$xi #get BW stutter setting
 xiFW = model$xiFW  #get FW stutter setting
 th0 <- mlefit$fit$thetahat #estimate params
 Sigma0 <- mlefit$fit$thetaSigma #estimated covariance of param estimtates
 np <- length(th0) #number of unknown parameters
 varnames <- names(mlefit$fit$thetahat) #variable names
 if(np!=ncol(Sigma0)) stop("Length of th0 and dimension of Sigma was not the same!")
 c <- mlefit$prepareC #returned from prepareC
 nC = model$nC #number of contributors
 usedeg <- !is.null(model$kit) #check for degradation
 
 #Prepare fixed params:
 nM = c$nM #number of markers to evaluate
 ATv = model$threshT
 pCv = model$prC
 lambdav = model$lambda
 fstv = model$fst
 
 #PREPEARING THE LIKELIHOOD OPTIMZATION (what parameters are provided?):
 useParamOther = rep(TRUE,3)  #index of parameters used (set NA if fixed)
 if(!usedeg) useParamOther[1] = FALSE #degrad not used
 if(!is.null(xi)) useParamOther[2] = FALSE #BW stutter not used
 if(!is.null(xiFW)) useParamOther[3] = FALSE  #FW stutter not used
 indexParamOther = rep(NA,3)
 if(any(useParamOther)) indexParamOther[useParamOther] = 1:sum(useParamOther) #init indices
 
 #The likelihood function taking theta-parameters (non-transformed params)
 logliktheta <- function(theta) {  
   if(any(theta<0)) return(-Inf) #Never consider negative params
   mixprop = as.numeric() #Length zero for 1 contributor
   if(nC>1) {
     mixprop = theta[1:(nC-1)] #extract parms
     if(any(mixprop>1)) return(-Inf) #Never consider mix prop above 1
   }
   
   muv = rep(theta[nC],nM)   #common param for each locus
   sigmav = rep(theta[nC+1],nM)  #common param for each locus
   
   #Prepare the remaining variables
   beta1 = 1.0 #set default values
   xiB = xi #set default values
   xiF = xiFW #set default values
   tmp = theta[ (nC+1) + indexParamOther ] #obtain params
   if(useParamOther[1]) {
     beta1 = tmp[1]
     if(beta1>1) return(-Inf)
   } 
   if(useParamOther[2]) {
     xiB = tmp[2]
     if(xiB>maxxi) return(-Inf)
   }
   if(useParamOther[3]) {
     xiF = tmp[3]
     if(xiF>maxxiFW) return(-Inf)
   } 
   betav = rep(beta1,nM) #common param for each locus
   xiBv = rep(xiB,nM) #common param for each locus
   xiFv = rep(xiF,nM) 
   
   #notice mixture proportions are unrestricted
   loglik = .C("loglikgammaC",as.numeric(0),c$nC,c$NOK,c$knownGind,as.numeric(mixprop),as.numeric(muv),as.numeric(sigmav),as.numeric(betav),as.numeric(xiBv),as.numeric(xiFv),as.numeric(ATv),as.numeric(pCv),as.numeric(lambdav),as.numeric(fstv),c$nReps,c$nM,c$nA,c$YvecLong,c$FvecLong,c$nTypedLong,c$maTypedLong,c$basepairLong,c$BWvecLong,c$FWvecLong,c$nPS,c$BWPvecLong,c$FWPvecLong,as.integer(maxThreads),as.integer(0),c$anyRel,c$relGind,c$ibdLong,PACKAGE="euroformix")[[1]]
   if(is.null(xi))  loglik <- loglik + log(model$pXi(xiB)) #weight with prior of xi
   if(is.null(xiFW))  loglik <- loglik + log(model$pXiFW(xiF)) #weight with prior of xiFW

   if(verbose) { #only show progressbar if verbose
     progcount <<- progcount + 1
     setTxtProgressBar(progbar,progcount)
   } 
   return(loglik) #weight with prior of tau and stutter.
 }
   
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
   postth[1,] <- th0 #init with proper value X[1,] #add noise into start point
   postlogL[1] <- loglik0 #init with proper value logliktheta(postth[1,]) #get start-likelihood  
   U <- runif(niter) #random numbers
   m <- 2 #counter for samples
   nacc <- 0  
   while(m<=niter) {
     postth[m,] <-  X[m,] + postth[m-1,] #proposed theta
     postlogL[m] <- logliktheta(theta=postth[m,])
     pr <- exp(postlogL[m]- postlogL[m-1]) #acceptance rate
     if(U[m]>pr) { #if not accepted, i.e. random prob too large (above pr)
      postth[m,] <-  postth[m-1,]
      postlogL[m] <- postlogL[m-1 ]
     } else {
      nacc <- nacc + 1
     }
     m <- m + 1 #update counter
   } #end while not done
  accrat <- nacc/(niter-1) #acceptance ratio (dont count first)
#  logpX <- logdmvnorm(postth,mean=th0,cholC=chol(Sigma0)) #insert with Normal-approx of post-th
  logpX <- logdmvnorm(t(postth),mean=th0,cholC=C) #insert with Normal-approx of post-th
  #plot(postlogL,ty="l")
  #plot(logpX,ty="l")
  
  logVals = logpX - postlogL #this is logged values intended to be "exped"
  offset = max(logVals) #find max value
  #margL <- 1/mean(exp(logVals - offset)) #estimated marginal likelihood
  logMargL <- log(length(logVals)) - offset - log(sum(exp(logVals-offset))) #estimated marginal likelihood (logged)

 if( verbose ) cat("\n") #new line
# nU <- nC-ret$nK #number of unknowns
# if(nU>1) { #if more than 1 unknown 
#  margL <- factorial(nU)*margL #get correct ML adjusting for symmetry
# } #end method
 colnames(postth) <- varnames  #save variable names
 return(list(logmargL=logMargL,posttheta=postth,postlogL=postlogL,logpX=logpX,accrat=accrat,MLE=mlefit$fit$thetahat,Sigma=mlefit$fit$thetaSigma,seed=seed))
} #end function

