#' @title calcTippet
#' @author Oyvind Bleka
#' @description Performing tippet/non-contributor analysis for a given set of hypotheses
#' @param tipIdx Index in condOrder which are replaced with a random man
#' @param mlefitHp A fitted object returned from calcMLE (under Hp)
#' @param mlefitHd A fitted object returned from calcMLE (under Hd)
#' @param niter Number of drawn non-contributors/iterations 
#' @param type type of Tippet analysis to conduct (MLE or INT)
#' @param LRobs The user can send an observed LR to superimpose to the plot
#' @param MLEopt MLE options (Not used since it uses what is stored in mlefit object)
#' @param INTopt INT options (must include reltol, maxEval, dev)
#' @param seed Seed used for reproducibility
#' @param verbose Whether printing the progress. Default is TRUE
#' @return returned non-contributor LR values and statistics
#' @export 

#type="MLE"; LRobs=NULL; INTopt=list(reltol=1, maxEval=1000,dev=2); seed = NULL; verbose=TRUE
calcTippet = function(tipIdx, mlefitHp, mlefitHd, niter=100, type="MLE", LRobs=NULL, MLEopt=NULL, INTopt=list(reltol=1, maxEval=1000,dev=2), seed = NULL, verbose=TRUE)  {
  if(!is.null(seed)) set.seed(seed)
  if(!type%in%c("MLE","INT")) stop("Type not supported!")
  
  #Obtain model variants
  if(is.null(LRobs) && type=="MLE") LRobs = (mlefitHp$fit$loglik-mlefitHd$fit$loglik)/log(10) #obtain observed LR for MLE (log10 scale)
  
  #Obtain settings:
  MLEopt = mlefitHd 
  maxThreads = mlefitHd$prepareC$maxThreads
  
  #Helpfunction to calculate the likelihood  
  calcLogLik = function(mod) { #noticed modified reference object
    #mod=mlefitHp$model
    #ALways calculate MLE based first      
    mle <- calcMLE(mod$nC,mod$samples,mod$popFreq,refDataRM,mod$condOrder,mod$knownRef, mod$kit, mod$DEG, mod$BWS, mod$FWS, 
                   mod$AT, mod$prC, mod$lambda, mod$fst, mod$knownRel, mod$ibd, mod$minF, mod$normalize, 
                   MLEopt$steptol, MLEopt$nDone, MLEopt$delta, MLEopt$difftol, MLEopt$seed, 
                   verbose=FALSE, priorBWS = mod$priorBWS, priorFWS = mod$priorFWS,
                   maxThreads=maxThreads, adjQbp=mod$adjQbp,resttol=MLEopt$resttol) 
    
    if(type=="MLE") {
      return(mle$fit$loglik)
      
      #Otherwise continue calculating 
    } else if(type=="INT") {
      lims = euroformix::getParamLimits(mle,dev=INTopt$dev)
      int <- calcINT(mle, lims$lower, lims$upper, reltol = INTopt$reltol, maxEval = INTopt$maxEval, scale = lims$scale, verbose=FALSE) 
      return(int$loglik)
    }
  }
  
  #Obtain Frequency information to Randomize non-contributors (provided under Hd)
  locs <- MLEopt$prepareC$markerNames #loci to evaluate
  fst = setNames(MLEopt$prepareC$fst,locs) #obtain fst per marker
  popFreq = MLEopt$model$popFreq
  refData = MLEopt$model$refData

  refKnownIdxHd = unique(c(which(mlefitHd$model$condOrder>0),mlefitHp$model$knownRef)) #index of known  individuals under Hd (contributors and known non-contributors)
  #refKnownIdxHp = setdiff(which(mlefitHp$model$condOrder>0),tipIdx) #index of known individuals under Hp (exlude tippet index=POI)
  refRelIdx = mlefitHd$model$knownRel #index of related individual
  ibd = mlefitHd$model$ibd
  calcUnderHd = any(fst>0) #check if re-calculation under Hd is needed
  if(!calcUnderHd && !is.null(ibd) && ibd[1]<1 && length(refRelIdx)==1 && refRelIdx==tipIdx) calcUnderHd = TRUE #must recalculate under Hd if POI is the related individual
  logLikHd = mlefitHd$fit$loglik #extract likelihood value under Hd
  
  # CALCULATE OUTCOME LIST OF RANDOM MAN (RM), with corresponding probabilities
  Glist <- list() 
  for(loc in locs) {
    refDataLoc = refData[[loc]] #keep a copy (assume correct structure already)
    if(!loc%in%names(refData)) refDataLoc = lapply(refData, function(x) x[[loc]]) #other structure
    freqAll = popFreq[[loc]] #get frequencies
    
    #Need to include missing alleles (avoid bug from v4.0.8):
    alleles = c(unlist(refDataLoc[refKnownIdxHd]),unlist(refDataLoc[refRelIdx])) #get alleles
    newA = alleles[!alleles%in%names(freqAll)] #new alleles
    if(length(newA)>0) {
      newA = setNames(rep(MLEopt$model$minF,length(newA)),newA)
      freqAll = c(freqAll,newA)
    }
    if(MLEopt$model$normalize) freqAll = freqAll/sum(freqAll)
    Glist[[loc]] = euroformix::calcGjoint(freq=freqAll,nU=1,fst=fst[loc],refK=unlist(refDataLoc[refKnownIdxHd]),refR=unlist(refDataLoc[refRelIdx]),ibd=ibd)
  } 
  #NOTICE THE KNOWN REFERENCE(S) GIVEN AS UNDER HD  
  if(verbose) print(paste0("Calculating ",niter," LR values...")) 
  
  #Initiate PROGRESS BAR:
  progcount = 1  #counter
  progbar <- txtProgressBar(min = 0, max = niter, style = 3) #create progress bar
  logLikHpRM <- logLikHdRM <- LRRM <- NULL
  for(m in seq_len(niter)) { #for each random individual from the population
    
    refDataRM = refData #copy reference data to include genotype of random profile
    #Draw random profile
    for(loc in locs) {
      GenoLoc = Glist[[loc]]
      randIdx = sample( seq_along(GenoLoc$Gprob),1,prob=GenoLoc$Gprob) #get random index
      GenoSim = as.character(GenoLoc$G[randIdx,]) #Obtain random genotype
      
      #Insert to correct structure
      if(loc%in%names(refData)) {
        refDataRM[[loc]][[tipIdx]] = GenoSim #include
      } else {
        refDataRM[[tipIdx]][[loc]] = GenoSim #include
      }
    }

    #CALCULATE UNDER Hp (and possibly under Hd)    
    logLikHp = calcLogLik(mlefitHp$model) #calculate under Hp
    if(calcUnderHd) logLikHd = calcLogLik(mod=mlefitHd$model) #calculate under Hd    

    #Store:    
    logLikHpRM[m] = logLikHp
    logLikHdRM[m] = logLikHd
    LRRM[m] = (logLikHp-logLikHd)/log(10) #convert to log10

    #Plot to figure after each x iterations:
    if(verbose && m%%(niter/10)==0)  plotTippet(LRRM[1:m],LRobs,mtxt=paste0(type,"-based"))

    progcount <- progcount + 1
    setTxtProgressBar(progbar,progcount) #update progress bar
  } #for each tippet
  stat = plotTippet(LRRM,lr0=LRobs,returnStatsOnly=TRUE)  #finally calculate statistics
  ret = list(LRobs=LRobs, LRRM=LRRM, logLikHpRM=logLikHpRM, logLikHdRM=logLikHdRM, type=type, stat=stat)  
  return(ret)
}