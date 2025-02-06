#' @title contLikSearch
#' @author Oyvind Bleka
#' @description Searches through a set of model configurations which optimizes the likelihood function
#'
#' @param NOC #The vector of number of contributors to search
#' @param modelDegrad Whether degradation should be modeled (requires kit to be specified): can be vector (FALSE,TRUE)
#' @param modelBWstutt Whether backward stutter should be modeled: can be vector (FALSE,TRUE)
#' @param modelFWstutt Whether forward stutter should be modeled: can be vector (FALSE,TRUE)
#' @param samples A List with samples which for each samples has locus-list elements with list elements adata and hdata. 'adata' is a qualitative (allele) data vector and 'hdata' is a quantitative (peak heights) data vector.
#' @param popFreq A list of allele frequencies for a given population.
#' @param refData Reference objects with (2-size) allele-vector given in list element [[i]][[s]].
#' @param condOrder condOrder under Hd. See calcMLE function for more details.
#' @param knownRefPOI Specify index of known POI under Hp (can only be one number)
#' @param prC A numeric for allele drop-in probability. Default is 0.
#' @param nDone Number of optimizations required providing equivalent results (same logLik value obtained)
#' @param threshT The detection threshold given. Used when considering probability of allele drop-outs.
#' @param fst The co-ancestry coefficient. Default is 0.
#' @param lambda Parameter in modeled peak height shifted exponential model. Default is 0.
#' @param delta Scaling of variation of normal distribution when drawing random startpoints. Default is 1.
#' @param kit Used to model degradation. Must be one of the shortnames of kit: check getKit()
#' @param verbose Whether printing optimization progress. Default is TRUE.
#' @param difftol Tolerance for being exact in log-likelihood value (relevant when nDone>1)
#' @param knownRel gives the index of the reference which the 1st unknown is related to.
#' @param ibd the identical by decent coefficients of the relationship (specifies the type of relationship)
#' @param minF The freq value included for new alleles (new alleles as potential stutters will have 0). Default NULL is using min.observed in popFreq.
#' @param normalize Whether normalization should be applied or not. Default is FALSE.
#' @param priorBWS Prior function for BWS-parameter. Flat prior on [0,1] is default.
#' @param priorFWS Prior function for FWS-parameter. Flat prior on [0,1] is default.
#' @param seed The user can set seed if wanted
#' @param maxThreads Maximum number of threads to be executed by the parallelization
#' @param alpha The significance level used for the envelope test (validMLEmodel). Default is 0.01
#' @param steptol Argument used in the nlm function for faster return from the optimization (tradeoff is lower accuracy).
#' @param adjQbp Indicate whether fragmenth length of Q-allele is based on averaged weighted with frequencies
#' @param resttol Restriction tolerance used to restrict genotype outcome. 0=No restriction, 1=Full restriction.
#' @return A list with fitted models, and table results for the defined model outcome
#' @export

contLikSearch = function(NOC=2:3, modelDegrad=FALSE,modelBWstutt=FALSE,modelFWstutt=FALSE, samples=NULL,popFreq=NULL,refData=NULL,condOrder=NULL,knownRefPOI=NULL, prC=0,nDone=3,threshT=50,fst=0,lambda=0,delta=1,kit=NULL,verbose=TRUE,difftol=0.01,knownRel=NULL,ibd=c(1,0,0),minF=NULL,normalize=TRUE,priorBWS=NULL,priorFWS=NULL,seed=NULL,maxThreads=0, alpha=0.01, steptol=1e-4, adjQbp = FALSE, resttol=1e-6){

  minNOC = 1 #minimum number of contributors to evaluate
  if(!is.null(condOrder)) minNOC = sum(condOrder>0) + 1 #minium NOC must exceed number of conditionals 
  if(minNOC>min(NOC)) {
    if(verbose) print("Minimum number of contributors exceeded the lowest assigned number of contributors.")
    NOC = NOC[minNOC<=NOC] 
    if(length(NOC)==0) return(NULL) #Return if nothing to search 
  }
  if(length(knownRefPOI)!=1) stop("Exactly one POI must be specified!")
  
  #Define conditional vectors:  
  condhp <- condhd <- condOrder  #copy conditional indices (assumed to be under Hd)
  condhp[knownRefPOI] = max(condhd) + 1 #add POI contributor as next conditonal contributor
  condhp[condhp>0] = 1:sum(condhp>0) #the indices will be sorted again (re-organized to follow same order as doing sole calculations)

  #Define known non-contributors
  knownNonContrHp = which(condhp==0)
  if(length(knownNonContrHp)==0) knownNonContrHp = NULL
  knownNonContrHd = c(knownNonContrHp,knownRefPOI)
  
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
  dontAllow = modoutcome[,3]==FALSE & modoutcome[,4]==TRUE #exclude this combination: FWS but not BWS
  modoutcome = modoutcome[!dontAllow,,drop=FALSE] #remove from outcome
  nMods = nrow(modoutcome) #number of models
  if(verbose) print(paste0("Evaluating ",nMods," model configurations. THIS MAY TAKE LONG TIME!"))

   booltxt = c("Y","N")
   fir = function(x)  ifelse(x,booltxt[1],booltxt[2]) #return Y/N (TRUE/FALSE)
   modtxt = apply(modoutcome,1,function(x) paste0("NOC=",x[1],"/DEG=",fir(x[2]),"/BWS=",fir(x[3]),"/FWS=",fir(x[4])))
   cn = c("logLik","adjLogLik","log10LR","MxPOI","SignifHp","SignifHd")
  
   hpfitList = list()
   hdfitList = list()
   outtable = matrix(nrow=nMods,ncol=length(cn)) #show table
   colnames(outtable) = cn
   rownames(outtable) = modtxt
   for(row in 1:nMods) { #traverse each models
    #row=1
    NOC = modoutcome[row,1] #get number of contributors to evaluate
    DEG = modoutcome[row,2]
    BWS = modoutcome[row,3]
    FWS = modoutcome[row,4]
  
    if(verbose) print(paste0("Evaluating model: NOC=",modoutcome[row,1]," DEG=",modoutcome[row,2]," BW=",modoutcome[row,3]," FW=",modoutcome[row,4]))
    if(FWS && !BWS) {
      if(verbose) print("Not possible to model Forward Stutter without assuming Backward Stutter!")
      next
    }
    
    if(verbose) print("Evaluating Hp:")
    hpfit <- hpfitList[[row]] <- calcMLE(NOC,samples,popFreq,refData,condhp,knownNonContrHp,kit,DEG,BWS,FWS,threshT,prC,lambda,fst,knownRel,ibd,minF,normalize,steptol,nDone,delta,difftol,seed,verbose,priorBWS,priorFWS,maxThreads,adjQbp,resttol)
    if(verbose) print("Evaluating Hd:")
    hdfit <- hdfitList[[row]] <- calcMLE(NOC,samples,popFreq,refData,condhd,knownNonContrHd,kit,DEG,BWS,FWS,threshT,prC,lambda,fst,knownRel,ibd,minF,normalize,steptol,nDone,delta,difftol,seed,verbose,priorBWS,priorFWS,maxThreads,adjQbp,resttol)
  
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
    POIcompontent = condhp[knownRefPOI] #obtain component of POI
    MxPOI = signif(hpfit$fit$thetahat2[POIcompontent],2) #extract estimated mix proportion
    
    newrow = c(round(c(loglik,adjloglik,LR),2),MxPOI,hpSignif,hdSignif) #new row in table
    outtable[row,] = newrow  #store data 
    if(verbose) print(outtable[1:row,,drop=FALSE]) #print table
  }

  return( list(modoutcome=modoutcome,hpfitList=hpfitList,hdfitList=hdfitList,outtable=outtable) )
}
