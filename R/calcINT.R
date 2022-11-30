
#' @title calcINT
#' @author Oyvind Bleka
#' @description Marginalizing the likelihood through numerical integration.
#' @details The procedure does numerical integration to approximate the marginal probability over the model parameters.
#' 
#' @param nC Number of contributors in model.
#' @param samples A List with samples which for each samples has locus-list elements with list elements adata and hdata. 'adata' is a qualitative (allele) data vector and 'hdata' is a quantitative (peak heights) data vector.
#' @param popFreq A list of allele frequencies for a given population.
#' @param lower Lower bounds of parameters. Must be in following order: mx1,..,mx_(nC-1),mu,sigma,beta,xi.
#' @param upper Upper bounds of parameters. Must be in following order: mx1,..,mx_(nC-1),mu,sigma,beta,xi.
#' @param refData Reference objects has locus-list element [[i]] with a list element 'r' which contains a 2 long vector with alleles for each references.
#' @param condOrder Specify conditioning references from refData (must be consistent order). For instance condOrder=(0,2,1,0) means that we restrict the model such that Ref2 and Ref3 are respectively conditioned as 2. contributor and 1. contributor in the model. 
#' @param knownRef Specify known non-contributing references from refData (index). For instance knownRef=(1,2) means that reference 1 and 2 is known non-contributor in the hypothesis. This affectes coancestry correction.
#' @param kit shortname of kit: Obtained from getKit()
#' @param DEG Boolean of whether Degradation model should be used
#' @param BWS Boolean of whether back-stutter model should be used
#' @param FWS Boolean of whether for-stutter model should be used
#' @param AT The analytical threshold given. Used when considering probability of allele drop-outs.
#' @param pC A numeric for allele drop-in probability. Default is 0.
#' @param lambda Parameter in modeled peak height shifted exponential model. Default is 0.
#' @param fst is the coancestry coeffecient. Default is 0.
#' @param knownRel gives the index of the reference which the 1st unknown is related to.
#' @param ibd the identical by decent coefficients of the relationship (specifies the type of relationship)
#' @param minF The freq value included for new alleles (new alleles as potential stutters will have 0). Default NULL is using min.observed in popFreq.
#' @param normalize Whether normalization should be applied or not. Default is FALSE.
#' @param priorBWS Prior function for xi-parameter (stutter). Flat prior on [0,1] is default.
#' @param priorFWS Prior function for xiFW-parameter (FW stutter). Flat prior on [0,1] is default.
#' @param reltol Required relative tolerance error of evaluations in integration routine. Default is 0.001.
#' @param scale used to make integrale calculateable for small numbers. For scale!=0, integrale must be scaled afterwards with exp(-scale) to be correct.
#' @param maxEval Maximum number of evaluations in the adaptIntegrate function. Default is 0 which gives an infinite limit.
#' @param verbose Whether printing limits to integrate over. Printing progress if maxEval>0. Default is TRUE.
#' @param maxThreads Maximum number of threads to be executed by the parallelization
#' @param adjQbp Indicate whether fragmenth length of Q-allele is based on averaged weighted with frequencies
#' @return ret A list(margL,deviation,nEvals,scale) where margL is Marginalized likelihood for hypothesis (model) given observed evidence, deviation is the confidence-interval of margL, nEvals is number of evaluations.
#' @export 
#' @references Hahn,T. (2005). CUBA - a library for multidimensional numerical integration. Computer Physics Communications, 168(2),78-95.
#' @keywords Marginalized likelihood


calcINT = function(nC,samples,popFreq, lower=NULL, upper=NULL, refData=NULL, condOrder = NULL, knownRef = NULL, kit=NULL,DEG=TRUE,BWS=TRUE,FWS=TRUE,
                   AT=50,pC=0.05,lambda=0.01,fst=0,knownRel=NULL,ibd=NULL,minF=NULL,normalize=TRUE, priorBWS=NULL, priorFWS=NULL, 
                   reltol=0.001, scale=0,maxEval=0,verbose=FALSE, maxThreads=0, adjQbp=FALSE) {
  
 if(is.null(maxEval)) maxEval <- 0
 if(is.null(lower) || is.null(upper)) stop("Not implemented for missing limits (do be done)")
 if(length(lower)!=length(upper)) stop("Length of integral limits differs")
 if(nC>1) {
   lower[1:(nC-1)] <- 0
   upper[1:(nC-1)] <- 1
 }

 start_time <- Sys.time()
 c = prepareC(nC,samples,popFreq, refData, condOrder, knownRef, kit,DEG,BWS,FWS,AT,pC,lambda,fst,knownRel,ibd,minF,normalize,adjQbp)
 #Other variables: basepairinfo, number of typed allelse  
 if(DEG && !c$useDEG) {
   print("Degradation had to be turned off since kitinfo was not found for selected kit.")
   DEG = FALSE
 }
 modTypes = c(DEG,BWS,FWS) #obtain model types
 if(length(lower)!= (nC+sum(modTypes)+1) ) stop("The length of the integral limits did not correpond with number of parameters!")
 if(DEG) upper[nC+2] = min(upper[nC+2],1) #restrict degradationif necessary
 if(BWS) {
   upper[nC + 2 + sum(DEG)] = min(upper[nC + 2 + sum(DEG)],1) #restrict BWS proportion if necessary
   if(FWS) upper[nC + 3 + sum(DEG)] = min(upper[nC + 3 + sum(DEG)],1) #restrict FWS proportion if necessary
 }
 
 #PREPARE FEEDING DATA
 if(verbose) print("Carrying out preparations for integration...")
 mod = Rcpp::Module( "mod",PACKAGE="euroformix" ) #load module
 obj = methods::new(mod$ExposedClass) #create object of class
 
 #Step 1: insert data to exposed class (filldata)
 obj$filldata(c$nStutterModels,c$nMarkers,c$nRepMarkers,c$nAlleles,c$startIndMarker_nAlleles,c$startIndMarker_nAllelesReps,c$peaks,c$freqs,c$dropinWeight, c$nTyped, c$maTyped, c$basepair,
              c$BWfrom, c$FWfrom, c$BWto, c$FWto, c$nPotStutters, c$startIndMarker_nAllelesTot, c$QalleleIndex, c$dropinProb, c$fst, c$AT, c$NOK, c$knownGind, c$relGind, c$ibd, as.integer(maxThreads)) 
 
 #Step 2: Indexing large matrix (doIndex)
 prepTime=system.time({
   obj$prepare(as.integer(nC))
 })[3]
 #sum((nA+nPS)*(nA*(nA+1)/2)^NOC) #number of evaluations
 if(verbose) print(paste0("Prep. done and took ",round(prepTime),"s. Start integration..."))
 
 #Inner likelihood function to integrate out
 liktheta <- function(theta) {
  param = .convBack(theta,nC, modTypes,isPhi=FALSE) #need to convert params back
  if(any(param<0)) return(0)
  loglik = obj$loglik(as.numeric(param))  #calculate likelihood
  if(!is.null(priorBWS) && param[nC+4]>0)  loglik <- loglik + log(priorBWS(param[nC+4])) #weight with prior of xi
  if(!is.null(priorFWS) && param[nC+5]>0)  loglik <- loglik + log(priorFWS(param[nC+5])) #weight with prior of xiFW
  likval <- exp(loglik+scale) #note the scaling given as parameter "scale".
  if(is.nan(likval)) return(0) #avoid NAN values
  
  if(verbose && maxEval>0) { #only show progressbar if verbose
    progcount <<- progcount + 1
    setTxtProgressBar(progbar,progcount)
  } 
  
  return(likval) #weight with prior of tau and stutter.
 }
 
 #DERIVED RESTRICTION FOR MIXTURE PROPORTIONS:
 nK = sum(condOrder>0) #number of conditionals
 nU <- nC-nK #number of unknowns
 if(nC==2 && nU==2) {
  lower[1] <- max(1/2,lower[1]) #restrict to 1/2-size
 }
 if(nC==3 && nU==3) { #restrict to 1/6-size
  lower[1] <-  max(1/3,lower[1])
  upper[2] <- min(1/2,upper[2])
 }
 if(nC==4 && nU==4) { #restrict to 1/12-size
  lower[1] <- max(1/4,lower[1])
  upper[2] <- min(1/2,upper[2])
  upper[3] <- min(1/3,upper[3])
 }
 if(nC==4 && nU==3) { #restrict to 1/2-size
  upper[3] <- min(1/2,upper[3])
 }
 if(nC==5 && nU==5) { #restrict to 1/20-size
  lower[1] <- max(1/5,lower[1])
  upper[2] <- min(1/2,upper[2])
  upper[3] <- min(1/3,upper[3])
  upper[4] <- min(1/4,upper[4])
 }
 if(nC==5 && nU==4) { #restrict to 1/3-size
  upper[3] <- min(1/2,upper[3])
  upper[4] <- min(1/3,upper[4])
 }
 if(nC==5 && nU==3) { #restrict to 1/2-size
  upper[4] <- min(1/2,upper[4])
 }
 if(nC==6 && nU==6) { #restrict to 1/25-size
  lower[1] <- max(1/5,lower[1])
  upper[2] <- min(1/2,upper[2])
  upper[3] <- min(1/3,upper[3])
  upper[4] <- min(1/4,upper[4])
  upper[5] <- min(1/5,upper[5])
 }
 if(nC==6 && nU==5) { #restrict to 1/4-size
  upper[3] <- min(1/2,upper[3])
  upper[4] <- min(1/3,upper[4])
  upper[5] <- min(1/4,upper[5])
 }
 if(nC==6 && nU==4) { #restrict to 1/3-size
  upper[4] <- min(1/2,upper[4])
  upper[5] <- min(1/3,upper[5])
 }
 if(nC==6 && nU==3) { #restrict to 1/2-size
  upper[5] <- min(1/2,upper[5])
 }

 #Get number of combinations which are used to scale the integral (cause of calculating symmetries):
 comb <- 1
 if(nC>1) {
   comb2 <- rep(1,nC-1) - (upper[1:(nC-1)]-lower[1:(nC-1)])
   comb <- round(1/prod(comb2[comb2>0]))
 }

#OVERVIEW OF RESTRICTIONS
#nU/nC 1  2  3  4  5  6
# 1	   1  1  1  1  1  1
# 2   NA  2  1  1  1  1
# 3   NA NA  6  2  2  2
# 4   NA NA NA 12  3  3
# 5   NA NA NA NA 20  4
# 6   NA NA NA NA NA 25

 if(verbose) {
   print(paste0("lower=",paste0(prettyNum(lower),collapse="/")))
   print(paste0("upper=",paste0(prettyNum(upper),collapse="/")))

   #Inititate progressbar if maxEval given
   progcount = 1  #counter
   if( maxEval>0 ) progbar <- txtProgressBar(min = 0, max = maxEval, style = 3) #create progress bar
 }
 
 
 foo <- cubature::adaptIntegrate(liktheta, lowerLimit = lower , upperLimit = upper , tol = reltol, maxEval=maxEval)#10000)
 if(verbose && maxEval>0) cat("\n") #skip line after progressbar
 obj$close() #free memory
 
 #Postvalues (need to scale values back)
 val <- foo$integral #estimated integral
 err <- foo$error #estimated relative error of integral
 
 #Take full integration into account (scales with comb)
 val <- comb*val
 err <- comb*err
 
 #Adjust with likelihood-scaling (done last)
 loglik = log(val)-scale #adjust to get loglik
 dev <- val + c(-1,1)*foo$error #unadjusted error
 logdev = log(dev) - scale #adjusted error (log-scale)
 dev = exp(-scale)*dev #same as exp(loglikError) (non-logged)
 return(list(loglik=loglik, logdeviation=logdev, margL=val,deviation=dev,nEvals=foo[[3]],scale=scale))
}

