#' @title qualLikMLE
#' @author Oyvind Bleka
#' @description Optimizing the likelihood function based on qualitative model (LRmix)
#'
#' @details The R-package forensim is used to calculate the likelihood value for given data and parameters.
#'
#' @param nC Number of contributors in model. Must be a constant.
#' @param samples A List with samples which for each samples has locus-list elements with list elements adata and hdata. 'adata' is a qualitative (allele) data vector and 'hdata' is a quantitative (peak heights) data vector.
#' @param popFreq A list of allele frequencies for a given population.
#' @param refData Reference objects with (2-size) allele-vector given in list element [[i]][[s]].
#' @param condOrder Specify conditioning references from refData (must be consistent order). For instance condOrder=(0,2,1,0) means that we restrict the model such that Ref2 and Ref3 are respectively conditioned as 2. contributor and 1. contributor in the model. 
#' @param knownRef Specify known non-contributing references from refData (index). For instance knownRef=(1,2) means that reference 1 and 2 is known non-contributor in the hypothesis. This affectes coancestry correction.
#' @param prC A numeric for allele drop-in probability. Can be a vector (must contain the marker names). Default is 0.
#' @param fst The co-ancestry coefficient. Can be a vector (must contain the marker names). Default is 0.
#' @param prDcontr assumed known dropout parameter for all contributors, NA means to be optimized. Must be a nC long vector if given.
#' @param prDcommon vector indicating which contributors should share common drop-out parameter. Assign integers to contributors with common parameters. NA means not optimized.
#' @param steptol Argument used in the nlm function for faster return from the optimization (tradeoff is lower accuracy).
#' @param maxIter Maximum number of iterations for the optimization
#' @param prDv0 Start values for fitting the drop-out probabilities (will be spanned if multidimensional)
#' 
#' @return ret A list(fit,model,nDone,delta,seed,prepareC) where fit is Maximixed likelihood elements for given model.
#' @export
#' @examples
#' \dontrun{
#' kit = "ESX17"
#' popfn = paste(path.package("euroformix"),"FreqDatabases",paste0(kit,"_Norway.csv")
#'	,sep=.Platform$file.sep)
#' evidfn = paste(path.package("euroformix"),"examples",paste0(kit,"_3p.csv"),
#'	sep=.Platform$file.sep)
#' popFreq = freqImport(popfn)[[1]] #obtain list with population frequencies
#' samples = sample_tableToList(tableReader(evidfn))
#' dat = prepareData(samples,popFreq=popFreq) #obtain data to use for analysis
#' fit = qualLikMLE(nC=2,samples=dat$samples,popFreq=dat$popFreq)
#' fit = qualLikMLE(nC=3,samples=dat$samples,popFreq=dat$popFreq,prDcommon=c(1,1,2))
#' }

qualLikMLE = function(nC,samples,popFreq,refData=NULL,condOrder=NULL,knownRef=NULL,prC=0.05,fst=0,prDcontr=NULL,prDcommon=NULL,steptol=1e-6,maxIter=100, prDv0= c(0.1,0.35,0.7)) {

  #Input checks:
  if(nC <= 0) stop("Number of contributors must be at least 1.")
  if(!is.null(condOrder)) {
    conds = condOrder[condOrder>0] #obtain conditional indices
    if( any(duplicated(conds)) ) stop("Wrong conditonal argument. Make sure to provide unique conditional indices.")
    if( nC < sum(condOrder>0) ) stop("The conditional number of references cannot exceed the number of contributors.")
    if( !is.null(knownRef) && knownRef%in%condOrder) stop("The known non-contributor assigned in knownRef argument can't be a contributor in condOrder argument.")
  } 
  
  sampleNames = names(samples)
  locs = names(popFreq) #loci to evaluate
  nM = length(locs) #number of markers/loci
  
  #Check dropin prob
  if(length(prC)==1) {
    pCv = rep(prC,nM) #common dropin prob for all markers
  } else {
    pCv = setVecRightOrder(prC, locs) 
  }
  
  #Check theta correction
  if(length(fst)==1) {
    fstv = rep(fst,nM) #common theta correction for all markers
  } else {
    fstv = setVecRightOrder(fst,  locs)
  }
  names(pCv) <- names(fstv) <- locs #insert locus names (correct order)
  
  nUnknown = nC - sum(condOrder>0) #number of unknowns
      
  #DATA IS PREPARED FOR OPTIMIZATION (for given hypothesis)
  #Note: STRINGS ARE ENCODED AS integers
  Evidlist <- popFreqQ <- Reflist <- list()
  for(loc in locs) { #for each marker
    Reflist[[loc]] <- list()

    #prepare allele frequencies (names)    
    popFreqQ[[loc]] = popFreq[[loc]]
    names(popFreqQ[[loc]]) <- 1:length(popFreq[[loc]])
    popA = names(popFreq[[loc]]) #obtain population alleles
    
    for(ref in names(refData[[loc]])) { #for each replicate
      #loc=locs[1]
      #ref=names(refData[[loc]])[1]
      adata <- refData[[loc]][[ref]]#obtain ref data
      if(length(adata)==0) {
        adata=NULL #is empty
      } else {
        adata = match(adata,popA) #obtain index
      }
      Reflist[[loc]][[ref]] <-adata #reference to store (index of popA)
    }
  
    adataSamples = NULL #store alleles across all samples
    for(s in 1:length(sampleNames)) { #for each replicate
      sample = sampleNames[s]
      adata <- samples[[sample]][[loc]]$adata #obtain evid data
      if(length(adata)==0) {
        adata=0 #is empty
      } else {
        adata = match(adata,popA) #obtain index
      }
      if( s>1 ) adata = c(0,adata) #separate replicates with zero
      adataSamples = c(adataSamples,adata)
    }
    Evidlist[[loc]] <- adataSamples #evidence to store (index of popA)
  }
  
  #Model calculation (make individual-specific dropout parameter possible):
  negloglik <- function(pD,returnLogLikPerMarkers=FALSE) {
    
    #transform back drop-out parameter
    if(!is.null(prDcommon)) {
      pDeachContr <-  1/(1+exp(-pD[prDcommon]))  #obtain dropout param for each contr (use corresponding index)
    } else {
      pDeachContr <- rep( 1/(1+exp(-pD)) , nC) #obtain dropout param for each contr
    }

    #Force input of fixed drop-out param
    if(!is.null(prDcontr)) {
      insdropoutind = !is.na(prDcontr) #index of where to insert known dropout vals
      pDeachContr[insdropoutind] = prDcontr[insdropoutind]
    }

    likval <- rep(1,length(locs))
    names(likval) = locs
    for(loc in locs) {
      fst0 = fstv[loc] #obtain marker specific setting (if set)
      pC0 = pCv[loc] #obtain marker specific setting (if set)
      condRefs = unlist(Reflist[[loc]][which(condOrder>0)]) #known contr
      knownNonCond = unlist(Reflist[[loc]][knownRef]) #known non-contr
      likval[which(loc==locs)] <- forensim::likEvid( Evidlist[[loc]],T=condRefs,V=knownNonCond,x=nUnknown,theta=fst0, prDHet=pDeachContr, prDHom=pDeachContr^2, prC=pC0, freq=popFreqQ[[loc]])
    }
    if(returnLogLikPerMarkers) {
      return( log(likval) )
    }  else {
      return( -sum(log(likval)) )
    }
  }
  
  #check if need to optimize (only if NULL or any element is NA)
  if(!is.null(prDcontr) && !any(is.na(prDcontr))) {
    opt = list()
    opt$logliki = negloglik( log(prDcontr/(1-prDcontr)),TRUE) #obtain loglik per markers
    opt$loglik <- sum(opt$logliki) #maximum (maximum likelihood)
    opt$pDhat = prDcontr #provide the fixed values 
    return(opt) #return from function without need to optimize
  }


  #in situation of having defined which dropout prob should be common for contributors
  if(!is.null(prDcommon)) {
    prDparUnique = unique(na.omit(prDcommon)) #obtain unique integers for 
    prDdim = length(prDparUnique)
    if(!all(prDparUnique==1:prDdim)) stop("Indices in argument prDcommon was wrongly defined (it must be 1:nparam).")
    tmpList = list() #help variable to obtain span of outcome
    for(i in 1:prDdim) tmpList[[i]] = prDv0 #insert to list
    prDm0 = expand.grid(tmpList) #expand grid to obtain startpoint outcome

    #Provide proper startvalues
    negLogLik0 = apply( log(prDm0/(1-prDm0)) ,1,negloglik) #loop through each outcome 
    indMin = which.min(negLogLik0) #index to use
    
    pD0 = unlist(prDm0[indMin,]) #starpoint to use
    opt <- nlm(negloglik,log(pD0/(1-pD0)), iterlim=maxIter, steptol=steptol)
  
  } else {

    #Provide proper startvalues
    negLogLik0 = Vectorize(negloglik)( log(prDv0/(1-prDv0)) )        #find suitable optim startpoint
    pD0 = prDv0[which.min(negLogLik0)] #starpoint to use
    opt <- nlm(Vectorize(negloglik),log(pD0/(1-pD0)), iterlim=maxIter, steptol=steptol)
    
  }
  
  opt$logliki = negloglik(opt$estimate,TRUE) #obtain loglik per markers
  opt$loglik <- -opt$min #maximum (maximum likelihood)
  opt$pDhat <- pD <- 1/(1+exp(-opt$estimate)) #obtian dropout-estimates 

  #################################
  #Insert dropout per contributors#
  #################################
  
  #transform back drop-out parameter
  if(!is.null(prDcommon)) {
    pDeachContr <- pD[prDcommon]  #obtain dropout param for each contr (use corresponding index)
  } else {
    pDeachContr <- rep(pD , nC) #obtain dropout param for each contr
  }
  
  #Force input of fixed drop-out param
  if(!is.null(prDcontr)) {
    insdropoutind = !is.na(prDcontr) #index of where to insert known dropout vals
    pDeachContr[insdropoutind] = prDcontr[insdropoutind]
  }
  
  opt$pDhatContr = pDeachContr #store dropout per contributors
  opt$model <- list(nC=nC,samples=samples,popFreq=popFreq,refData=refData,condOrder=condOrder,knownRef=knownRef,prC=pCv,fst=fstv,prDcontr=prDcontr,prDcommon=prDcommon)
  opt$steptol=steptol
  opt$maxIter=maxIter
  opt$prDv0=prDv0
  
  return(opt)
} #end function
  