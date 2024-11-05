#' @title calcINT
#' @author Oyvind Bleka
#' @description Marginalizing the likelihood through numerical integration.
#' @param mlefit Fitted object using calcMLE
#' @param lower Lower bounds of parameters. Must be in following order: mx1,..,mx_(nC-1),mu,sigma,beta,xi.
#' @param upper Upper bounds of parameters. Must be in following order: mx1,..,mx_(nC-1),mu,sigma,beta,xi.
#' @param reltol Required relative tolerance error of evaluations in integration routine.
#' @param scale used to avoid underflow (should be the maximum likelihood value)
#' @param maxEval Maximum number of evaluations in the adaptIntegrate function. 
#' @param verbose Whether printing limits to integrate over. Printing progress if maxEval>0
#' @return ret A list(loglik,logdeviation,nEvals,scale) where loglik is Marginalized likelihood, logdeviation is the error-interval of loglik, nEvals is the total number of likelihood evaluations.
#' @export 
#' @references Hahn,T. (2005). CUBA - a library for multidimensional numerical integration. Computer Physics Communications, 168(2),78-95.
#' @keywords Marginalized likelihood


calcINT = function(mlefit, lower, upper, reltol=1,scale=0,maxEval=1000,verbose=FALSE) {
  if(is.null(lower) || is.null(upper)) stop("Not implemented for missing limits (do be done)")
  if(length(lower)!=length(upper)) stop("Length of integral limits differs")

  tol2integrate = 0.01 #tolerance for integrating a variable (width of integral)
  model <- mlefit$model
  priorBWS = model$priorBWS #obtain prior stutter model
  priorFWS = model$priorFWS #obtain prior stutter model
  modTypes <- mlefit$modelType #obtain model types
  c <- mlefit$prepareC #returned from prepareC
  nC = model$nC #number of contributors
  
  #PREPARE:
  if(verbose) print("Carrying out preparations for Integral...")
  obj = prepareCobj(c, verbose) #use wrapper function to obtain C++ pointer
  
  restTime=system.time({ #prepare restriction (can take a while)
    preSearched = .preSearch(obj,c,mlefit$thetaPresearch,mlefit$resttol,priorBWS, priorFWS)
  })[3]
  if(verbose) {
    print(paste0("Presearch done and took ",round(restTime),"s. Start integrate..."))
    print(paste0("Percentage of genotype outcome to be used: ",round(preSearched$restprop*100,1),"% (tol=",mlefit$resttol,")"))
  }
  nEvals = sum(sapply(preSearched$presearchFitList,function(x) nrow(x))) #store number of evaluations used
  
  #prepare restriction rule for Mx  
  nKinship = as.integer(c$hasKinship) #number of kinship specifications
  if(nKinship>0) stop("Cannot (yet) calculate INT based LR where kinship is specified!")
  nKnown = c$nC - c$nU + nKinship#number of conditionals
  nUnknown = c$nU - nKinship#number of Mx params to be restricted
  
  #The likelihood function taking theta-parameters (non-transformed params)
 # mxLower = pmax(0,lower[1:(nC-1)]) #avoid lower than zero
 # mxUpper = pmin(1,upper[1:(nC-1)]) #avoid higher than one
  other_lower = pmax(0,lower[nC:length(lower)]) #avoid lower than zero
  other_upper = upper[nC:length(upper)]
  #Also avoid that DEG,StuttProp params are higher than one
  #if(length(other_upper)>2) other_upper[-c(1:2)] = pmin(1,other_upper[-c(1:2)])
  
  mxLims = NULL
  if(nC>1) { #only relevant for having at least 2 contributors
    mxLims = cbind(0,1/(2:nC)) #obtain Mx limits in case of unknowns
    if(nKnown>0) { #if conditionals
      mxLimsKnown = t(replicate(nKnown,c(0,1))) #get limits for knowns
      mxLims = rbind(mxLimsKnown,mxLims)[1:(nC-1),,drop=FALSE]
    }
    mxLower = pmax(mxLims[,1],lower[1:(nC-1)])
    mxUpper = pmin(mxLims[,2],upper[1:(nC-1)])
    mxLims = cbind(mxLower,mxUpper) #construct limit
  }
  #Init an own integration function
  myInt = function(fun,L,U,...) {
    diff = U-L
    if(diff<tol2integrate) { #dont integrate if less than this (use middle)
      middle = (U+L)/2 #insert "middle" as known variable
      int = fun(middle,...)
    } else {
      int = cubature::adaptIntegrate(fun,L,U, tol=2,...)[[1]]
      #int = integrate(Vectorize(fun),L,U, rel.tol=1,...)[[1]]
    }
    return(int)
  }
  
  
  logLik <- function(par) {
    nEvals <<- nEvals + 1
    param = .convBack(par,nC, modTypes,isPhi=FALSE) #need to convert params back
    llv = obj$loglik(as.numeric(param) )
    if(!is.null(priorBWS) && param[nC+4]>0)  llv <- llv + log(priorBWS(param[nC+4])) #weight with prior of xi
    if(!is.null(priorFWS) && param[nC+5]>0)  llv <- llv + log(priorFWS(param[nC+5])) #weight with prior of xiFW
    return(llv)
  }

#DEFINING DIFFERENT INTEGRAL FUNCTION FOR DIFFERENT HYPOTHESIS
  LikFun = function(par) { #helpfunction to obtain (Scaled) likelihood value
    return(exp(logLik(par)-scale))
  }
  

  
####DONE DEFINING OUTCOME OF INTEGRAL
  
  if(verbose) {
    print(paste0("lower=",paste0(round(lower,2),collapse="/")))
    print(paste0("upper=",paste0(round(upper,2),collapse="/")))
  }
  
  #Indicate calculation situations
  calcIntMx = NULL
  if(nC==1) calcIntMx = .calcInt1p
  if(nC==2) {
    if(nKnown==0) calcIntMx = .calcInt2p0K
    if(nKnown>0) calcIntMx = .calcInt2p1K
  }
  if(nC==3) {
    if(nKnown==0) calcIntMx = .calcInt3p0K
    if(nKnown==1) calcIntMx = .calcInt3p1K
    if(nKnown>1) calcIntMx = .calcInt3p2K
  }
  if(nC==4) {
    if(nKnown==0) calcIntMx = .calcInt4p0K
    if(nKnown==1) calcIntMx = .calcInt4p1K
    if(nKnown==2) calcIntMx = .calcInt4p2K
    if(nKnown>2) calcIntMx = .calcInt4p3K
  }
  if(nC==5) {
    if(nKnown==0) calcIntMx = .calcInt5p0K
    if(nKnown==1) calcIntMx = .calcInt5p1K
    if(nKnown==2) calcIntMx = .calcInt5p2K
    if(nKnown==3) calcIntMx = .calcInt5p3K
    if(nKnown>3) calcIntMx = .calcInt5p4K
  }
  if(is.null(calcIntMx)) {
    print("WARNING:Cannot calculate integral for at least 6 contributors!")
    return(NULL)
  } 
  
 
  
  #Create a wrapper function to print out progress
  progbar = NULL #init progbar
  progcount = 1  #counter
  intfunWrapper = function(thetaOther,stuttParam=NULL,...) {
    if(!is.null(stuttParam)) thetaOther = c(thetaOther,stuttParam) #always append value
    if(progcount==1 && verbose) print("Start integrating...")
    timeMxInt = system.time({ #calculate time for one call
      val = calcIntMx(LikFun,myInt,mxLims,thetaOther,...)
    })[3]
    if(verbose && maxEval>0) { #only show progressbar if verbose
      if(progcount==1) {
        expectedTimeProgress = .secondToTimeformat(timeMxInt*maxEval)
        print(paste0("Upper time for integration: ", expectedTimeProgress, " (HH:MM:SS):")) 
        
        #Inititate progressbar
        progbar <<- txtProgressBar(min = 0, max = maxEval, style = 3) #create progress bar
      }
      setTxtProgressBar(progbar,progcount)
    } 
    progcount <<- progcount + 1
#    print(progcount)
    return(val)
  }
 
  ######################################## 
  #PERFORM INTEGRATION ('Outer integral')#
  #Special situation: Check if the integral width of the stutter parameters is too narrow
  
  # Identify stutter indices based on model types
  stuttIdx = NULL
  if(modTypes[2]) stuttIdx =  sum(modTypes[1]) + 3 #also count PHexp/PHvar variables
  if(modTypes[3]) stuttIdx <- c(stuttIdx, sum(modTypes[1:2]) + 3)

  # Determine if we need to skip integration
  outer_diff <- other_upper - other_lower  # Obtain width of integrals
  dontIntegrate <- outer_diff[stuttIdx] < tol2integrate
  
  if(any(dontIntegrate)) { #there was a special situation where we want to avoid integrate
    #in case of having both BW and FW model ..and to Integrate on FWstutt but not BWstutt (Weird!)
    if(length(stuttIdx)==2 && dontIntegrate[1] && !dontIntegrate[2]) { #
      stop("Progress with integration is stalled because the BW-stutter model was not informative, whereas the FW-stutter model was-an unexpected outcome.")
    }
    # Update stutter indices to skip integration for narrow stutter parameters
    stuttIdx <- stuttIdx[dontIntegrate]
    outer_MLE <- (other_upper + other_lower) / 2  # Obtain MLE
    #mlefit$fit$thetahat2
    integration_result <- cubature::adaptIntegrate(intfunWrapper, lowerLimit=other_lower[-stuttIdx] ,upperLimit=other_upper[-stuttIdx], tol=reltol, maxEval=maxEval, stuttParam=outer_MLE[stuttIdx])
  } else { #integrate as normal
    integration_result <- cubature::adaptIntegrate(intfunWrapper, lowerLimit=other_lower ,upperLimit=other_upper, tol=reltol, maxEval=maxEval)
  }
  if(verbose && maxEval>0)  cat("\n") #ensure skipping line after computations
  
  #intFun = cubature::hcubature
  #integration_result <- cubature::adaptIntegrate(intfunWrapper, lowerLimit=other_lower ,upperLimit=other_upper, tol=intFunRelTol, maxEval=maxEval)
  obj$close() #free memory
  
  #FINALIZING RESULTS
  # Retrieve estimated integral and error (scaled to avoid overflow)
  estimated_integral <- integration_result$integral
  estimated_error <- integration_result$error
  
  #Obtain the correct likelihood value (Adjust with scaling)
  log_likelihood <- log(estimated_integral) + scale  

  # Calculate unadjusted error range
  deviation_range <- estimated_integral + c(-1, 1) * estimated_error

  # Adjust the error range for log-scale likelihood (Adjust with scaling)
  log_deviation_range <- log(deviation_range) + scale  # Log-scale adjusted error range
 # exp(scale)*deviation_range (un-logged)
  
  return(list(loglik=log_likelihood, logdeviation=log_deviation_range, nEvals=nEvals,scale=scale))
}


