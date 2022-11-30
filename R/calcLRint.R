#' @title calcLRint
#' @author Oyvind Bleka
#' @description Calculating marginal likelihood based on integration
#' @details The procedure is a wrapper for the numerical integration function.
#' @param mlefitHp A fitted object returned from contLikMLE (under Hp)
#' @param mlefitHd A fitted object returned from contLikMLE (under Hd)
#' @param reltol Required relative tolerance error of evaluations in integration routine. Default is 0.001.
#' @param maxEval Maximum number of evaluations in the adaptIntegrate function. Default is 0 which gives an infinite limit.
#' @param dev Scaling relative to 2*SE (i.e. 97.5 percentile). Should cover most of post
#' @param verbose Whether printing simulation progress. Default is TRUE
#' @return returned inferred LR with error range (reflects relative error)
#' @export 

calcLRint = function(mlefitHp, mlefitHd, reltol=0.001,maxEval=0,dev=3,verbose=TRUE) {
  
  #helpfunction for integral
  calcfun = function(mle) {
    lims = euroformix::getParamLimits(mle,dev) #obtain values under Hd if not provided otherwise
    model = mle$model
    time = system.time({
      int = euroformix::calcINT(nC=model$nC,samples=model$samples,popFreq=model$popFreq, lower=lims$lower, upper=lims$upper, 
                    refData=model$refData, condOrder=model$condOrder, knownRef=model$knownRef, kit=model$kit,
                    DEG=model$DEG,BWS=model$BWS,FWS=model$FWS, AT=model$AT,pC=model$prC,lambda=model$lambda,
                    fst=model$fst,knownRel=model$knownRel,ibd=model$ibd,minF=model$minF,normalize=model$normalize, 
                    priorBWS=model$priorBWS, priorFWS=model$priorFWS, reltol=reltol, scale=lims$scale,maxEval=maxEval,verbose=verbose, 
                    maxThreads=mle$maxThreads, adjQbp=model$adjQbp)
    })[3]
    int$time = time
    int$paramLimits = lims
    if(verbose) {
      print(paste0("Time usage: ", .secondToTimeformat(time),"(HH:MM:SS)"))
      print(paste0("Number of evaluations: ",int$nEvals))
      print(paste0("log(Lik)=",int$loglik))
      print(paste0("logLik-error=[",paste0(int$logdev,collapse=","),"]"))
      #print(paste0("Lik=",exp(int$loglik))) #use small numbers here
    }
    return(int)
  }

  
  if(verbose) print("Calculations done for Hp...")
  intHp = calcfun(mlefitHp)
  
  
  if(verbose) print("Calculations done for Hd...")
  intHd = calcfun(mlefitHd)
  
  #Post-calculations
  lr = (intHp$loglik-intHd$loglik)
  #LRerror <- range(c(intHp$deviation/intHd$deviation,intHp$deviation/rev(intHd$deviation))) #get deviation interval of LR
  lrdev = c(intHp$logdev[1]-intHd$logdev[2], intHp$logdev[2]-intHd$logdev[1])
  
  ret = list(LR=exp(lr),log10LR=lr/log(10),LRerror=exp(lrdev),log10LRerror=lrdev/log(10),
             calcHp=intHp, calcHd=intHd, reltol=reltol,maxEval=maxEval,dev=dev)
  return(ret)
} #end function

