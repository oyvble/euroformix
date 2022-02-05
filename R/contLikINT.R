
#' @title contLikINT
#' @author Oyvind Bleka
#' @description contLikINT marginalizes the likelihood through numerical integration.
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
#' @param xi A numeric giving stutter-ratio if it is known. Default is NULL, meaning it is integrated out.
#' @param prC A numeric for allele drop-in probability. Default is 0.
#' @param reltol Required relative tolerance error of evaluations in integration routine. Default is 0.001.
#' @param threshT The detection threshold given. Used when considering probability of allele drop-outs.
#' @param fst is the coancestry coeffecient. Default is 0.
#' @param lambda Parameter in modeled peak height shifted exponential model. Default is 0.
#' @param pXi Prior function for xi-parameter (stutter). Flat prior on [0,1] is default.
#' @param kit shortname of kit: Obtained from getKit()
#' @param scale used to make integrale calculateable for small numbers. For scale!=0, integrale must be scaled afterwards with exp(-scale) to be correct.
#' @param maxEval Maximum number of evaluations in the adaptIntegrate function. Default is 0 which gives an infinite limit.
#' @param knownRel gives the index of the reference which the 1st unknown is related to.
#' @param ibd the identical by decent coefficients of the relationship (specifies the type of relationship)
#' @param xiFW A numeric giving FW stutter-ratio if it is known.Default is 0, meaning stutter is not used.
#' @param pXiFW Prior function for xiFW-parameter (FW stutter). Flat prior on [0,1] is default.
#' @param maxThreads Maximum number of threads to be executed by the parallelization
#' @param verbose Whether printing limits to integrate over. Printing progress if maxEval>0. Default is TRUE.
#' @return ret A list(margL,deviation,nEvals,scale) where margL is Marginalized likelihood for hypothesis (model) given observed evidence, deviation is the confidence-interval of margL, nEvals is number of evaluations.
#' @export 
#' @references Hahn,T. (2005). CUBA - a library for multidimensional numerical integration. Computer Physics Communications, 168(2),78-95.
#' @keywords Marginalized likelihood


contLikINT = function(nC,samples,popFreq,lower,upper,refData=NULL,condOrder=NULL,knownRef=NULL,xi=NULL,prC=0,reltol=0.01,threshT=50,fst=0,lambda=0,pXi=function(x)1,kit=NULL,scale=0,maxEval=0,knownRel=NULL,ibd=c(1,0,0),xiFW=0,pXiFW=function(x)1,maxThreads=32,verbose=TRUE){
 if(is.null(maxEval)) maxEval <- 0
 if(length(lower)!=length(upper)) stop("Length of integral limits differs")
 np = length(lower)

 c <- prepareC(nC,samples,popFreq,refData,condOrder,knownRef,kit,knownRel,ibd,fst,incS=is.null(xi) || xi>0,incFS=is.null(xiFW) || xiFW>0)
 usedeg  <- !is.null(kit) #boolean whether modeling degradation TRUE=YES, FALSE=NO
 
 #Prepare fixed params:
 nM = c$nM #number of markers to evaluate
 
 #Check AT:
 if(length(threshT)==1) {
   ATv = rep(threshT,nM) #common detection threshold for all markers
 } else {
   ATv = setVecRightOrder(threshT,  c$locs) #get right order of vector 
 }
 
 #Check dropin prob
 if(length(prC)==1) {
   pCv = rep(prC,nM) #common dropin prob for all markers
 } else {
   pCv = setVecRightOrder(prC,  c$locs) 
 }
 
 #Check hyperparam DI
 if(length(lambda)==1) {
   lambdav = rep(lambda,nM) #common hyperparam DI for all markers
 } else {
   lambdav = setVecRightOrder(lambda,  c$locs) 
 }
 
 #Check theta correction
 if(length(fst)==1) {
   fstv = rep(fst,nM) #common theta correction for all markers
 } else {
   fstv = setVecRightOrder(fst,  c$locs)
 }
 names(ATv) <- names(pCv) <- names(lambdav) <- names(fstv) <- c$locs #insert locus names (correct order)
 
 #PREPEARING THE LIKELIHOOD OPTIMZATION (what parameters are provided?):
 useParamOther = rep(TRUE,3)  #index of parameters used (set NA if fixed)
 if(!usedeg) useParamOther[1] = FALSE #degrad not used
 if(!is.null(xi)) useParamOther[2] = FALSE #BW stutter not used
 if(!is.null(xiFW)) useParamOther[3] = FALSE  #FW stutter not used
 indexParamOther = rep(NA,3)
 if(any(useParamOther)) indexParamOther[useParamOther] = 1:sum(useParamOther) #init indices
 
 liktheta <- function(theta) {
  mixprop = as.numeric() #Length zero for 1 contributor
  if(nC>1) mixprop = theta[1:(nC-1)] #extract parms
  muv = rep(theta[nC],nM)   #common param for each locus
  sigmav = rep(theta[nC+1],nM)  #common param for each locus

  #Prepare the remaining variables
  beta1 = 1.0 #set default values
  xiB = xi #set default values
  xiF = xiFW #set default values
  tmp = theta[ (nC+1) + indexParamOther ] #obtain params
  if(useParamOther[1]) beta1 = tmp[1]
  if(useParamOther[2]) xiB = tmp[2]
  if(useParamOther[3]) xiF = tmp[3]
  betav = rep(beta1,nM) #common param for each locus
  xiBv = rep(xiB,nM) #common param for each locus
  xiFv = rep(xiF,nM) 
    
  loglik = .C("loglikgammaC",as.numeric(0),c$nC,c$NOK,c$knownGind,as.numeric(mixprop),as.numeric(muv),as.numeric(sigmav),as.numeric(betav),as.numeric(xiBv),as.numeric(xiFv),as.numeric(ATv),as.numeric(pCv),as.numeric(lambdav),as.numeric(fstv),c$nReps,c$nM,c$nA,c$YvecLong,c$FvecLong,c$nTypedLong,c$maTypedLong,c$basepairLong,c$BWvecLong,c$FWvecLong,c$nPS,c$BWPvecLong,c$FWPvecLong,as.integer(maxThreads),as.integer(0),c$anyRel,c$relGind,c$ibdLong,PACKAGE="euroformix")[[1]]
  if(is.null(xi))  loglik <- loglik + log(pXi(xiB)) #weight with prior of xi
  if(is.null(xiFW))  loglik <- loglik + log(pXiFW(xiF)) #weight with prior of xiFW
  likval <- exp(loglik+scale) #note the scaling given as parameter "scale".
  
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
 
# NOC1  2  3  4  5  6
# 1  1  1  1  1  1  1
# 2 NA  2  1  1  1  1
# 3 NA NA  6  2  2  2
# 4 NA NA NA 12  3  3
# 5 NA NA NA NA 20  4
# 6 NA NA NA NA NA 25

 if(verbose) {
   print(paste0("lower=",paste0(prettyNum(lower),collapse="/")))
   print(paste0("upper=",paste0(prettyNum(upper),collapse="/")))

   #Inititate progressbar if maxEval given
   progcount = 1  #counter
   if( maxEval>0 ) progbar <- txtProgressBar(min = 0, max = maxEval, style = 3) #create progress bar
 }
 
 
 foo <- cubature::adaptIntegrate(liktheta, lowerLimit = lower , upperLimit = upper , tol = reltol, maxEval=maxEval)#10000)
 val <- foo$integral
 dev <- val + c(-1,1)*foo$error
 nEvals <- foo[[3]]
 val <- comb*val
 dev <- comb*dev
 return(list(margL=val,deviation=dev,nEvals=nEvals,scale=scale))
}

