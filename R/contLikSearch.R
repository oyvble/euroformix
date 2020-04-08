#' @title contLikSearch
#' @author Oyvind Bleka
#' @description contLikSearch search through the models  optimizes the likelihood function of the DNA mixture model 
#' @details Function calls the likelihood function implemented in C++ which uses the package Boost and paralellisation using OpenMP
#'
#' @param NOC #The vector of number of contributors to search
#' @param modelDegrad Boolean of whether degradation should be modeled (requires kit to be specified): can be vector (FALSE,TRUE)
#' @param modelBWstutt Boolean of whether backward stutter should be modeled: can be vector (FALSE,TRUE)
#' @param modelFWstutt Boolean of whether forward stutter should be modeled: can be vector (FALSE,TRUE)
#' @param samples A List with samples which for each samples has locus-list elements with list elements adata and hdata. 'adata' is a qualitative (allele) data vector and 'hdata' is a quantitative (peak heights) data vector.
#' @param popFreq A list of allele frequencies for a given population.
#' @param refData Reference objects with (2-size) allele-vector given in list element [[i]][[s]].
#' @param condOrder Specify conditioning references from refData (must be consistent order). For instance condOrder=(0,2,1,0) means that we restrict the model such that Ref2 and Ref3 are respectively conditioned as 2. contributor and 1. contributor in the model. 
#' @param knownRefPOI Specify known POI under Hp, under non-contributing references (under Hd) from refData (index). For instance knownRef=(1,2) means that reference 1 and 2 is known non-contributor in the hypothesis. This affectes coancestry correction.
#' @param prC A numeric for allele drop-in probability. Default is 0.
#' @param nDone Number of optimizations required providing equivalent results (same logLik value obtained)
#' @param threshT The detection threshold given. Used when considering probability of allele drop-outs.
#' @param fst The co-ancestry coefficient. Default is 0.
#' @param lambda Parameter in modeled peak height shifted exponential model. Default is 0.
#' @param pXi Prior function for xi-parameter (stutter). Flat prior on [0,1] is default.
#' @param delta Scaling of variation of normal distribution when drawing random startpoints. Default is 1.
#' @param kit Used to model degradation. Must be one of the shortnames of kit: check getKit()
#' @param verbose Boolean whether printing optimization progress. Default is TRUE.
#' @param maxIter Maximum number of iterations for the optimization. Default is 100.
#' @param knownRel gives the index of the reference which the 1st unknown is related to.
#' @param ibd the identical by decent coefficients of the relationship (specifies the type of relationship)
#' @param pXiFW Prior function for xiFW-parameter (FW stutter). Flat prior on [0,1] is default.
#' @param seed The user can set seed if wanted
#' @param maxThreads Maximum number of threads to be executed by the parallelization
#' @param alpha The significance level used for the envelope test (validMLEmodel). Default is 0.01
#' @return ret A list(fit,model,nDone,delta,seed,prepareC) where fit is Maximixed likelihood elements for given model.
#' @export
contLikSearch = function(NOC=2:3, modelDegrad=FALSE,modelBWstutt=FALSE,modelFWstutt=FALSE, samples=NULL,popFreq=NULL,refData=NULL,condOrder=NULL,knownRefPOI=NULL, prC=0,nDone=2,threshT=50,fst=0,lambda=0,pXi=function(x)1,delta=1,kit=NULL,verbose=TRUE,maxIter=100,knownRel=NULL,ibd=c(1,0,0),pXiFW=function(x)1,seed=NULL,maxThreads=32, alpha=0.01){

  minNOC = 1 #minimum number of contributors to evaluate
  if(!is.null(condOrder)) minNOC = sum(condOrder>0) + 1 #minium NOC must exceed number of conditionals 
  if(minNOC>min(NOC)) {
    if(verbose) print("Minimum number of contributors exceeded the lowest assigned number of contributors.")
    NOC = NOC[minNOC<=NOC] 
    if(length(NOC)==0) return(NULL) #Return if nothing to search 
  }
  
  #fix conditional vectors:  
  condhp <- condhd <- condOrder  #copy conditional indices
  condhp[knownRefPOI] = max(condhd) + 1 #add POI contributor
  
  #Spanning model outcome (based on model settings):
  modList = list()
  modList[[1]] = NOC #minNOC:nC #range of number of contributors
  
  if(any(modelDegrad) && is.null(kit)) {
    print("Kit not selected, degradation model will not be considered")
    modelDegrad = FALSE #degradation turned off
  }
  modList[[2]] <- modelDegrad
  modList[[3]] <- modelBWstutt
  modList[[4]] <- modelFWstutt
  
  #SPAN MODEL OUTCOME
  modoutcome = expand.grid(modList) #get model outcome
  colnames(modoutcome) = c("NOC","DEG","BWstutt","FWstutt")
  modoutcome = modoutcome[order(modoutcome[,1]),,drop=FALSE] #reorder
  nMods = nrow(modoutcome) #number of models
  if(verbose) print(paste0("Evaluating ",nMods," model configurations. THIS MAY TAKE LONG TIME!"))

 booltxt = c("Y","N")
 fir = function(x)  ifelse(x,booltxt[1],booltxt[2]) #return Y/N (TRUE/FALSE)
 modtxt = apply(modoutcome,1,function(x) paste0("NOC=",x[1],"/DEG=",fir(x[2]),"/BWS=",fir(x[3]),"/FWS=",fir(x[4])))
 cn = c("logLik","adjLogLik","log10LR","MxPOI","SignifHp","SignifHd")
 outtable = matrix(nrow=nMods,ncol=length(cn)) #show table
 colnames(outtable) = cn
 rownames(outtable) = modtxt
 
 hpfitList = list()
 hdfitList = list()
  for(row in 1:nMods) { #traverse each models
    NOC = modoutcome[row,1] #get number of contributors to evaluate
    kit0 = NULL #default if not modeled
    if(modoutcome[row,2]) kit0 = kit #get kitname  #DEGRADATION
    xi <- xiFW <- 0 #default if not modeled 
    if(modoutcome[row,3]) xi = NULL #BACKWARD STUTTER
    if(modoutcome[row,4]) xiFW = NULL #FORWARD STUTTER
    
    if(verbose) print(paste0("Evaluating model: NOC=",modoutcome[row,1]," DEG=",modoutcome[row,2]," BW=",modoutcome[row,3]," FW=",modoutcome[row,4]))
    if(verbose) print("Evaluating Hp:")
    hpfit <- hpfitList[[row]] <- contLikMLE(NOC,samples,popFreq,refData,condhp,NULL,       xi,prC,nDone,threshT,fst,lambda,pXi,delta,kit0,verbose,maxIter,knownRel,ibd,xiFW,pXiFW,seed,maxThreads)
    if(verbose) print("Evaluating Hd:")
    hdfit <- hdfitList[[row]] <- contLikMLE(NOC,samples,popFreq,refData,condhd,knownRefPOI,xi,prC,nDone,threshT,fst,lambda,pXi,delta,kit0,verbose,maxIter,knownRel,ibd,xiFW,pXiFW,seed,maxThreads)

    hpSignif <- hdSignif <- NA #default is no values
    if(!is.infinite(hpfit$fit$loglik)) {
      hpvalid = validMLEmodel(hpfit,kit,alpha=alpha,createplot=FALSE,verbose=FALSE) #check model
      hpSignif = sum(hpvalid$Significant)
    }
    if(!is.infinite(hdfit$fit$loglik)) {
      hdvalid = validMLEmodel(hdfit,kit,alpha=alpha,createplot=FALSE,verbose=FALSE) #check model
      hdSignif = sum(hdvalid$Significant)
    }
    LR = (hpfit$fit$loglik - hdfit$fit$loglik)/log(10) #obtain log10 LR
    loglik = hdfit$fit$loglik
    adjloglik = loglik - length(hdfit$fit$thetahat) #penalty by number of params
    MxPOI = signif(hpfit$fit$thetahat2[knownRefPOI],2) #extract estimated mix proportion
    
    newrow = c(round(c(loglik,adjloglik,LR),2),MxPOI,hpSignif,hdSignif) #new row in table
    outtable[row,] = newrow  #store data 
    if(verbose) print(outtable[1:row,,drop=FALSE]) #print table
  }

  return( list(modoutcome=modoutcome,hpfitList=hpfitList,hdfitList=hdfitList,outtable=outtable) )
}
