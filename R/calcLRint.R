#' @title calcLRint
#' @author Oyvind Bleka
#' @description Calculating marginal likelihood based on integration
#' @details The procedure is a wrapper for the numerical integration function for obtaining Bayes Factor
#' @param mlefitHp A fitted object returned from calcMLE (under Hp)
#' @param mlefitHd A fitted object returned from calcMLE (under Hd)
#' @param reltol Required relative tolerance error of evaluations in integration routine
#' @param maxEval Maximum number of evaluations in the integration routine function. 
#' @param dev Scaling relative to 2*SE (i.e. 97.5 percentile). Should cover most of posterior distribution
#' @param verbose Whether printing simulation progress. 
#' @return returned inferred LR with error range (reflects relative error)
#' @export 

calcLRint = function(mlefitHp, mlefitHd, reltol=1,maxEval=1000,dev=2,verbose=TRUE) {

  calcfun = function(mlefit) {
    lims =  getParamLimits(mlefit,dev) #obtain values under Hd if not provided otherwise
    time = system.time({
      int = calcINT(mlefit, lower=lims$lower, upper=lims$upper,reltol=reltol,scale=lims$scale,maxEval=maxEval,verbose=verbose)
    })[3]
    int$time = time
    int$paramLimits = lims
    if(verbose) {
      print(paste0("Time usage: ", .secondToTimeformat(time),"(HH:MM:SS)"))
      print(paste0("Number of evaluations: ",int$nEvals))
      print(paste0("log(Lik)=",round(int$loglik,2)))
      print(paste0("logLik-error=[",paste0(round(int$logdev,2),collapse=","),"]"))
      #print(paste0("Lik=",exp(int$loglik))) #use small numbers here
    }
    return(int)
  }

  
  if(verbose) print("Performing calculations for Hp...")
  intHp = calcfun(mlefitHp)
  
  
  if(verbose) print("Performing calculations for Hd...")
  intHd = calcfun(mlefitHd)
  
  #Post-calculations
  lr = (intHp$loglik-intHd$loglik)
  #LRerror <- range(c(intHp$deviation/intHd$deviation,intHp$deviation/rev(intHd$deviation))) #get deviation interval of LR
  lrdev = c(intHp$logdev[1]-intHd$logdev[2], intHp$logdev[2]-intHd$logdev[1])
  nEvals = c(intHp$nEvals,intHd$nEvals) #obtain number of evaluations
     
  ret = list(LR=exp(lr),log10LR=lr/log(10),LRerror=exp(lrdev),log10LRerror=lrdev/log(10),
             calcHp=intHp, calcHd=intHd, reltol=reltol,nEvals=nEvals,dev=dev)
  return(ret)
} #end function

