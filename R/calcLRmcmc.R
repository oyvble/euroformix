#' @title calcLRmcmc
#' @author Oyvind Bleka
#' @description A function for calculating LR from MCMC simulations
#' @details Returns a conservative (quantile) LR or a Full Bayesian based LR
#' @param mlefitHp A fitted object returned from contLikMLE (under Hp)
#' @param mlefitHd A fitted object returned from contLikMLE (under Hd)
#' @param niter Number of samples in the MCMC-sampling.
#' @param delta A numerical parameter to scale with the covariance function Sigma. Default is 2. Will be tuned automtically.
#' @param quantile The quantile used to report conservative LR 
#' @param seed The user can set seed if wanted
#' @param accRate aimed acceptance rate
#' @param accRateTol Difference tolerance regarding acceptance rate
#' @param niterTune Number of samples for tuning delta in the MCMC-sampling.
#' @param diffSeed Seed difference
#' @param verbose Whether progress should be printed
#' @param traceplot Whether to add a traceplot with 95\% CI
#' @param mcmcObjList A list with hp/hd earlier run with contLikMCMC output
#' @return A list with different LR statistics from the MCMC based simulations
#' @export

#niter=5000;delta=2;quantile=0.05;seed=NULL; accRate=0.25;accRateTol=0.1; niterTune=200; diffSeed=999; verbose=TRUE; traceplot=TRUE;mcmcObjList=NULL
calcLRmcmc = function(mlefitHp, mlefitHd,niter=2000,delta=2,quantile=0.05,seed=NULL, accRate=0.25,accRateTol=0.1, niterTune=200, diffSeed=999, verbose=TRUE, traceplot=TRUE, mcmcObjList=NULL) {
  #PERFORM CALIBRATING OF DELTA BEFORE RUNNING ALL SAMPLE
  #Tweak delta to find  suitable acceptance rate:
  LRmle = exp(mlefitHp$fit$loglik - mlefitHd$fit$loglik) #obtain MLE based LR
  deltaInput = delta #store input delta
  
  if(is.null(mcmcObjList)) { #only calibrate if first time
    if(verbose) print("Calibrating MCMC simulator...")
    while(TRUE) {  #only run if no earlier results
      if(verbose) print(paste0("Check with delta=",delta))
      hpmcmc <- contLikMCMC(mlefitHp,niter=niterTune,delta=delta)
      acc0 = hpmcmc$accrat #obtain acceptance rate
      if(verbose) print(paste0("Acceptance rate=",acc0))
      if( abs(acc0-accRate)<accRateTol) break
      
      if(acc0==0) {
        scaleAcc = 0.5 #reduce sampling variation by 1/2 if none is accepted
      } else {
        scaleAcc = acc0/accRate #obtain scale between accepted and Golden
      }
      delta = delta*scaleAcc #update delta
    }
    if(verbose) print(paste0("Tuned delta=",delta))
  }
  gfun = function(x) quantile(x,quantile) #this is function to obtain result from

  if(verbose) print("Sampling under Hp...")
  hpmcmc <- contLikMCMC(mlefitHp,niter=niter,delta=delta,seed=seed,mcmcObj=mcmcObjList$hp)
  
  if(verbose) {
    print(paste0("Estimated integral: logLik=",hpmcmc$logmargL))
    print("Sampling under Hd...")
  }
  seed2 = seed + diffSeed #same seed diff as in EFM
  if(length(seed2)==0) seed2 = NULL
  hdmcmc <- contLikMCMC(mlefitHd,niter=niter,delta=delta,seed=seed2,mcmcObj=mcmcObjList$hd) 
  if(verbose) print(paste0("Estimated integral: logLik=",hdmcmc$logmargL))
  delta = hdmcmc$delta #be sure that latest delta is the one tuned.
  
  #Post-evaluation:
  log10LRdistr <- (hpmcmc$postlogL - hdmcmc$postlogL)/log(10) #calculate log10LR
  LRcons  <- gfun(log10LRdistr)
  niter = length(log10LRdistr) #get updated length

  #Estimation of bayesian based LR:
  LRbayes <- exp(hpmcmc$logmargL-hdmcmc$logmargL)#/log(10) #convert log to log10
  if(is.infinite(LRbayes)) LRbayes <- NaN #invalid number due to zero Hd
  
  #estimate number of effective samples
  #Perform bootstrap to estimate uncertainty estimate of quantile:
  bootCI = function(x,nrep=1000,adjusted=TRUE,alpha=0.05) {
    size = length(x)
    if(adjusted) {
      acfobj = acf(x,lag.max=100,plot=FALSE)
      size = round( length(x)/(1+2*sum(acfobj$acf[-1])) ) #efficient sample saize
    }    
    qboot = replicate(nrep, gfun(sample(x,size,replace=TRUE))) 
    return(quantile(qboot, c(alpha/2, 1-alpha/2)) )
  }
  qCI = bootCI(log10LRdistr,nrep=10000) #final CI

  #Obtain trace-plot:
  if(traceplot) {
    if(verbose) print("Creating traceplot...")
    
    #obtain trace of Bayes Factor
    margLogLik = function(mcmc) {
      logVals = mcmc$logpX - mcmc$postlogL #this is logged values intended to be "exped"
      offset = mcmc$offset
      return( log(seq_along(logVals)) - offset - log(cumsum(exp(logVals-offset)))) #estimated marginal likelihood (logged)
    }
    BayesFactorTrace =  (margLogLik(hpmcmc) - margLogLik(hdmcmc))/log(10)
    #plot(1:length(BayesFactorTrace),BayesFactorTrace,ty="l")
    
    niterTrace1 = floor(seq(min(niter,500),niter,l=100)) #used to show smoothed trace
    niterTrace2 = floor(seq(min(niter,1000),niter,l=10)) #used to show 95% CI (less places)
    LRconsTrace  <- sapply(niterTrace1, function(x) gfun(log10LRdistr[1:x]))
    LRconsCITrace  <- sapply(niterTrace2, function(x) bootCI(log10LRdistr[1:x]))
    ylim = range(c(LRconsCITrace,BayesFactorTrace,log10(LRmle)+1))
    #LRconsCITrace= LRconsTrace
    plot(niterTrace1,LRconsTrace,ty="l",ylim=ylim,ylab="log10LR",xlab="Number of iterations",main="Trace plot (MCMC)",col=4)
    mtext(paste0("Quantile for conservative: ",quantile))
    
    #Explanation:  
    #txt2 = paste0("95% CI=[",paste0(round(qCI,2),collapse=","),"]")
    legtxt = c("MLE based","BayesFactor (Trace)","BayesFactor (Final)","Conservative (Trace)","Conservative (Final)","Conservative (95% CI)")
    legend("topleft",legtxt,col=c(1,2,2,4,4,4),lty=c(2,1,2,1,2,1),cex=0.8,pch=c(rep(NA,5),19))
    abline(h=log10(LRmle),lty=2,col=1)
    abline(h=qCI,lty=3,col=4)
    abline(h=LRcons,lty=2,col=4)
    for(i in 1:2) lines(niterTrace2,LRconsCITrace[i,],ty="o",col=4,pch=19)
    lines(seq_along(BayesFactorTrace),BayesFactorTrace,col=2) 
    if(!is.nan(LRbayes)) abline(h=log10(LRbayes),col=2,lty=2)
  }
  

  mcmcList = list(LRbayes=LRbayes, log10LRcons=LRcons, log10LRbayes=log10(LRbayes), 
                  LRcons=10^LRcons, LRmle=LRmle, log10LRdistr=log10LRdistr,
                  niter=niter,delta=deltaInput,quantile=quantile,seed=seed,
                  accRate=accRate,accRateTol=accRateTol,deltaTuned=delta,
                  niterTune=niterTune, diffSeed=diffSeed, log10LRconsCI = qCI,
                  mcmcObjList = list(hp=hpmcmc,hd=hdmcmc))
 return(mcmcList)
}