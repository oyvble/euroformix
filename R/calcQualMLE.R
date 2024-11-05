#' @title calcQualMLE
#' @author Oyvind Bleka
#' @description Optimizing the likelihood function based on qualitative model (LRmix).
#' @param nC Number of contributors in model. Must be a constant.
#' @param samples A List with samples which for each samples has locus-list elements with list elements adata and hdata. 'adata' is a qualitative (allele) data vector and 'hdata' is a quantitative (peak heights) data vector.
#' @param popFreq A list of allele frequencies for a given population.
#' @param refData Reference objects with (2-size) allele-vector given in list element [[locus]][[sample]].
#' @param condOrder Specify conditioning references from refData (must be consistent order). For instance condOrder=(0,2,1,0) means that we restrict the model such that Ref2 and Ref3 are respectively conditioned as 2. contributor and 1. contributor in the model. 
#' @param knownRef Specify known non-contributing references from refData (index). For instance knownRef=(1,2) means that reference 1 and 2 is known non-contributor in the hypothesis. This affectes coancestry correction.
#' @param prC A numeric for allele drop-in probability. Can be a vector (must contain the marker names). Default is 0.
#' @param fst The co-ancestry coefficient. Can be a vector (must contain the marker names). Default is 0.
#' @param prDcontr assumed known dropout parameter for all contributors, NA means to be optimized. Must be a nC long vector if given.
#' @param prDcommon vector indicating which contributors should share common drop-out parameter. Assign integers to contributors with common parameters. NA means not optimized.
#' @param steptol Argument used in the nlm function for faster return from the optimization (tradeoff is lower accuracy).
#' @param prDv0 Start values for fitting the drop-out probabilities (will be spanned if multidimensional)
#' @param maxIter Maximum number of iterations for nlm
#' @return ret A list(fit,model,nDone,delta,seed,prepareC) where fit is Maximixed likelihood elements for given model.
#' @export
#' @examples
#' \dontrun{
#' AT = 50 #analytical threshold
#' sep0 = .Platform$file.sep
#' popfn = paste(path.package("euroformix"),"FreqDatabases",paste0(kit,"_Norway.csv"),sep=sep0)
#' evidfn = paste(path.package("euroformix"),"examples",paste0(kit,"_3p.csv"),sep=sep0)
#' reffn = paste(path.package("euroformix"),"examples",paste0(kit,"_refs.csv"),sep=sep0)
#' popFreq = freqImport(popfn)[[1]] #population frequencies
#' samples = sample_tableToList(tableReader(evidfn)) #evidence samples
#' refData = sample_tableToList(tableReader(reffn)) #reference sample
#' dat = prepareData(samples,refData,popFreq,threshT=AT) #needed for qual method
#' condOrder = c(1,2,0) #assuming C1=ref1,C2=ref2
#' logLik1 = calcQualMLE(2,dat$samples,dat$popFreq,dat$refData,condOrder)$loglik
#' logLik2 = calcQualMLE(3,dat$samples,dat$popFreq,dat$refData,condOrder, prDcommon=c(1,1,2))$loglik
#' }

calcQualMLE = function(nC,samples,popFreq,refData=NULL,condOrder=NULL,knownRef=NULL,prC=0.05,fst=0,prDcontr=NULL,prDcommon=NULL,steptol=1e-6, prDv0= c(0.1,0.35,0.7), maxIter = 10000) {

  #Input checks:
  if(nC <= 0) stop("Number of contributors must be at least 1.")
  if(!is.null(condOrder)) {
    condRefsIdx = which(condOrder>0) #get references to condition on
    conds = condOrder[condRefsIdx] #obtain conditional indices
    if( any(duplicated(conds)) ) stop("Wrong conditonal argument. Make sure to provide unique conditional indices.")
    if( nC < sum(condOrder>0) ) stop("The conditional number of references cannot exceed the number of contributors.")
    if( !is.null(knownRef) && length(knownRef)>0 && length(condRefsIdx)>0 && knownRef%in%condRefsIdx) stop("The known non-contributor assigned in knownRef argument can't be a contributor in condOrder argument.")
  } 
  
  sampleNames = names(samples)
  locs_samples = unique(unlist(lapply(samples,names)))
  locs_pop = names(popFreq) #loci to evaluate
  locs = intersect(locs_samples,locs_pop) #use those in common only
  nLocs = length(locs) #number of markers/loci

  #Check dropin prob
  if(length(prC)==1) {
    pCv = setNames(rep(prC,nLocs),locs) #common dropin prob for all markers
  } else {
    pCv = setVecRightOrder(prC, locs) 
  }
  
  #Check theta correction
  if(length(fst)==1) {
    fstv = setNames(rep(fst,nLocs),locs) #common theta correction for all markers
  } else {
    fstv = setVecRightOrder(fst,  locs)
  }
  names(pCv) <- names(fstv) <- locs #insert locus names (correct order)
  
  nCond = sum(condOrder>0) #number of conditionals
  nUnknown = setNames(rep(NA,length(locs)),locs) #create vector for number of unknowns
  
  #DATA IS PREPARED FOR OPTIMIZATION (for given hypothesis)
  #Note: STRINGS ARE ENCODED AS integers
  Evidlist <- popFreqQ <- Reflist <- list()
  for(loc in locs) { #for each marker
  #  loc=locs[6]
    Reflist[[loc]] <- list()
    
    #prepare allele frequencies (names)    
    popFreqQ[[loc]] = popFreq[[loc]]
    names(popFreqQ[[loc]]) <- 1:length(popFreq[[loc]])
    popA = names(popFreq[[loc]]) #obtain population alleles (BEFORE ENCODING!)
    refNames = names(refData[[loc]])
    
    nCondRefs = 0
    for(r in seq_along(refNames)) { #for each reference
      ref = refNames[r]
      #loc=locs[1]
      #ref=names(refData[[loc]])[1]
      adata <- unlist(refData[[loc]][[ref]]) #obtain ref data
      if(length(adata)==0) {
        adata=as.character() #is empty
      } else {
        if(condOrder[r]>0) nCondRefs = nCondRefs + 1 #
        adata = match(adata,popA) #obtain index
      }
      Reflist[[loc]][[ref]] <-adata #reference to store (index of popA)
    }
    
    #Obtain number of non-empty references in condOrder 
    nUnknown[loc] = nC - nCondRefs
      
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
      condRefs <- knownNonCond <- NULL
      if(nCond>0) condRefs = unlist(Reflist[[loc]][condRefsIdx]) #known contr
      if(length(knownRef)>0) knownNonCond = unlist(Reflist[[loc]][knownRef]) #known non-contr
      likval[which(loc==locs)] <- calcQual( Evidlist[[loc]], condRefs, knownNonCond, nUnknown[loc], fst0, pDeachContr, pC0, popFreqQ[[loc]])
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
  opt$prDv0=prDv0
  
  return(opt)
} #end function
  