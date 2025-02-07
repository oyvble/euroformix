#' @title calcMLE
#' @author Oyvind Bleka
#' @description Optimizes the likelihood function (MLE)
#' @param nC Number of contributors in model.
#' @param samples A list of samples (evidence) with structure [[samplename]][[locus]] = list(adata,...)
#' @param popFreq A list of allele frequencies for a given population.
#' @param refData Reference objects has locus-list element with a list element 'r' which contains a 2 long vector with alleles for each references.
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
#' @param steptol Argument used in the nlm function for faster return from the optimization (tradeoff is lower accuracy).
#' @param nDone Number of optimizations required providing equivalent results (same logLik value obtained)
#' @param delta NOT IN USE Scaling of variation of normal distribution when drawing random startpoints. Default is 1.
#' @param difftol NOT IN USE Tolerance for being exact in log-likelihood value (relevant when nDone>1)
#' @param seed NOT IN USE The user can set seed if wanted 
#' @param verbose Whether to print progress of optimization
#' @param priorBWS Prior function for BWS-parameter. Flat prior on [0,1] is default.
#' @param priorFWS Prior function for FWS-parameter. Flat prior on [0,1] is default.
#' @param maxThreads Maximum number of threads to be executed by the parallelization
#' @param adjQbp Indicate whether fragmenth length of Q-allele is based on averaged weighted with frequencies
#' @param resttol Restriction tolerance used to restrict genotype outcome. 0=No restriction, 1=Full restriction.
#' @return Fitted maximum likelihood object 
#' @export 
#' @examples
#' \dontrun{
#' kit = "ESX17"
#' AT = 50 #analytical threshold
#' sep0 = .Platform$file.sep
#' popfn = paste(path.package("euroformix"),"FreqDatabases",paste0(kit,"_Norway.csv"),sep=sep0)
#' evidfn = paste(path.package("euroformix"),"examples",paste0(kit,"_3p.csv"),sep=sep0)
#' reffn = paste(path.package("euroformix"),"examples",paste0(kit,"_refs.csv"),sep=sep0)
#' popFreq = freqImport(popfn)[[1]] #population frequencies
#' samples = sample_tableToList(tableReader(evidfn)) #evidence samples
#' refData = sample_tableToList(tableReader(reffn)) #reference sample
#' plotEPG2(samples,kit,refData,AT)
#' condOrder = c(1,2,0) #assuming C1=ref1,C2=ref2
#' mlefit = contLikMLE(3,samples,popFreq,refData,condOrder,kit=kit)
#' plotTopEPG2(mlefit)
#' }


calcMLE = function(nC,samples,popFreq, refData=NULL, condOrder = NULL, knownRef = NULL, kit=NULL,DEG=TRUE,BWS=TRUE,FWS=TRUE,
  AT=50,pC=0.05,lambda=0.01,fst=0,knownRel=NULL,ibd=NULL,minF=NULL,normalize=TRUE,steptol = 1e-4,nDone=3,delta=1,difftol=0.01,
  seed=NULL,verbose=FALSE,priorBWS=NULL,priorFWS=NULL, maxThreads=0, adjQbp = FALSE, resttol=1e-6) {

  start_time <- Sys.time()
  if(!is.null(seed)) set.seed(seed) #set seed if provided

  #Preparing variables to fill into C++ structure
  c = prepareC(nC,samples,popFreq, refData, condOrder, knownRef, kit,BWS,FWS,AT,pC,lambda,fst,knownRel,ibd,minF,normalize, adjQbp)
  if(DEG && !c$useDEG) {
    print("Degradation had to be turned off since kitinfo was not found for selected kit.")
    DEG = FALSE
  }
  modTypes = setNames(c(DEG,BWS,FWS),c("DEG","BWS","FWS")) #obtain model types
  c$modTypes = modTypes #also store in c object
  c$maxThreads = as.integer(maxThreads) #stored only here
    
  #Prefitting data based on the model for sum of the peak heights  to find proper startvalues for MLE 
  sumY <- meanbp <- rep(NA,c$nMarkers)
  startind1 = 0 #start index (no reps)
  startindR = 0 #start index (all reps)
  for(i in seq_len(c$nMarkers)) { #for each marker
    nAreps = c$nAlleles[i]*c$nRepMarkers[i] #number of alleles over all reps
    ind1 = startind1 + c(1:c$nAlleles[i]) #obtain index of PH vector to use
    indR = startindR + c(1:nAreps) #obtain index of PH vector to use
    sumY[i] <- sum(c$peaks[indR])/c$nRepMarkers[i] #take sum of the peak heights (average over number of reps)
    meanbp[i] <- mean(c$basepair[ind1]) #take sum of the peak heights
    startind1 =  tail(ind1,1)
    startindR =  tail(indR,1)
  }
 # plot(meanbp,sumY)
  ncond = sum(condOrder>0) #numbers to conidtion on
  c$nU = nC-ncond #obtain number of unknowns (also store in c object)
  np = (nC+1) + DEG + BWS + FWS  #obtain number of parameters
  if(!DEG) meanbp = NULL
  #obtain expected startpoints of (PHexp,PHvar,DegSlope)
  th0 <- fitgammamodel(y=sumY,x=meanbp,niter=10,delta=1,offset=0,scale=1) 
  thetaPresearch = th0
  thetaPresearch[2] = max(thetaPresearch[2],0.2) #use minimum 0.2 in PHvar in presearch
  if(!DEG) thetaPresearch = c(thetaPresearch,1)
  
  #PREPARE FEEDING DATA INTO C++ structure
  if(verbose) print("Carrying out preparations for optimizing...")
  obj = prepareCobj(c, verbose) #use wrapper function to obtain C++ pointer
 
#  nTotEvals = (c$nAlleles+c$nPotStutters)*(c$nAlleles*(c$nAlleles+1)/2)^c$nC #number of evaluations

#NEW FROM v4.1.0 #############################################################################    
#The program will perform a pre-search of likely paramters and store the highest obtained genotype likelihoods
#The idea is to restrict the genotype outcome to only use the most important ones (a tolerance is used to define this)

  restTime=system.time({ #prepare restriction (can take a while)
    preSearched = .preSearch(obj,c,thetaPresearch,resttol,priorBWS, priorFWS)
  })[3]
  if(verbose) print(paste0("Presearch done and took ",round(restTime),"s. Start optimizing..."))
  presearchFitList = preSearched$presearchFitList
  restprop = preSearched$restprop
  nEvals = sum(sapply(presearchFitList,function(x) nrow(x))) #store number of evaluations used
  if(verbose) print(paste0("Percentage of genotype outcome to be used: ",round(restprop*100,1),"% (threshold=",resttol,")"))
  
 #END NEW#############################################################################    
  
  #function for calling on C-function: Must convert "real domain" values back to model params 
  negloglikYphi = function(phi,progressbar=TRUE) { #assumed order: mixprop(1:C-1),mu,sigma,beta,xi
    param = .convBack(phi,nC, modTypes) #must be full parameter vector
    loglik = obj$loglik(as.numeric(param) ) 
    paramStutters = param[length(param) + c(-1,0)] #obtain stutter params (always last two)
    if(any(paramStutters>=1)) return(Inf) #dont allow above 1
    if(!is.null(priorBWS) && paramStutters[1]>0) loglik = loglik + log( priorBWS(paramStutters[1]) )
    if(!is.null(priorFWS) && paramStutters[2]>0) loglik = loglik + log( priorFWS(paramStutters[2]) )
    
    if(progressbar) {
      progcount <<- progcount + 1 #update counter
      if(verbose) setTxtProgressBar(progbar,progcount)  #only show progressbar if verbose(and if decided to show)
    }
    #print(loglik)
    return(-loglik) #weight with prior of stutter.
  }  

  #Narrow down candidate searches after preSearch
  #Obtain %-wise best solutions in pre-search (per stutter types)
  truncatePercentageLevel = 0.01 #keep all with at least 1%
  #  table =  presearchFitList[[1]]
  bestPreSearchParamsList = lapply(presearchFitList, function(table) {
    maxLikInds = 1
    if(nrow(table)>1) {
      loglik0vec = table[,ncol(table)]
      lik0vec = exp(loglik0vec-max(loglik0vec)) #avoid overflow by subtract max
      lik0vec = lik0vec/sum(lik0vec)
      maxLikInds = which(lik0vec>truncatePercentageLevel) #all with more than 1% likelihood
    }
    #barplot(lik0vec)
    return(table[maxLikInds,,drop=FALSE])#,-ncol(table)])
  })
  
  #Arrange all search sets into one matrix and remove duplicated Mx-sets:
  #The top "nDone" sets will be kept (zero-stutter values are lowest on list unless they provide unique Mx sets)
  bestPreSearchParams = do.call("rbind",bestPreSearchParamsList)
  isDup = duplicated(bestPreSearchParams[,1:nC,drop=FALSE])
  bestPreSearchParams = bestPreSearchParams[!isDup,,drop=FALSE]
  if(nrow(bestPreSearchParams)>1) {
    ord = order(bestPreSearchParams[,ncol(bestPreSearchParams)],decreasing=TRUE)
    nbestUse = max(1,min(nDone,length(ord))) #here nDone is used to steer maximal number of optimizations
    bestPreSearchParams = bestPreSearchParams[ord[1:nbestUse],,drop=FALSE]
  }
  
  #Evaluate one time to get timer:
  maxIterProgress = (np-2)*100 #maximum iterations in progressbar
  logLiks = bestPreSearchParams[,ncol(bestPreSearchParams)]
  maxIdx = which.max(logLiks)
  maxL <- logLiks[maxIdx]
  if(length(maxL)==0) maxL = -Inf #insert in case of none found
  thetaStart = .getThetaUnknowns(bestPreSearchParams[maxIdx,],nC,modTypes) #obtain start value of theta
  
  phi0 = .getPhi(thetaStart,nC,modTypes) #Note:Remove restrictive
  timeOneCall = system.time({ #estimate the time for calling the likelihood function one time 
    negLikVal <- negloglikYphi(phi=phi0,FALSE)   #check if start value was accepted
  })[3] #obtain time in seconds
  valdiff = abs(negLikVal+maxL) #get loglik diff
  if(!is.nan(valdiff) && !is.na(valdiff) && valdiff>difftol) print("WARNING: The performed restriction may give inaccurate likelihood values. You should reduce the restriction threshold.")
  
  #Show upper expected time:
  expectedTimeProgress0 = timeOneCall*maxIterProgress
  expectedTimeProgress = .secondToTimeformat(expectedTimeProgress0)
  showProgressBar = expectedTimeProgress0 > 10 #show if more than 10s
  if(verbose && showProgressBar) print(paste0("Expected (upper) time per optimization is ", expectedTimeProgress, " (HH:MM:SS):")) 
  
  #TRAVERSE
  maxPhi <- rep(NA,np) #Set as NULL to recognize if able to be estimated
  for(startParamIdx in seq_len(nrow(bestPreSearchParams)))  {
#    startParamIdx=1
    thetaStart = .getThetaUnknowns(bestPreSearchParams[startParamIdx,],nC,modTypes) #obtain params to use

    if(verbose) {
      print(paste0("Performing optimization ",startParamIdx,"/",nrow(bestPreSearchParams),"...."))
      #print(paste0("Start value used: ",paste0(round(thetaStart,2),collapse="/")))
    }

    #PREPARE PROGRESSBAR
    progcount  = 1 #progress counter for optimization (reset for each sucessful optimization)
    if(verbose && showProgressBar) {
      progbar <- txtProgressBar(min = 0, max = maxIterProgress, style = 3) #create progress bar)
    }

    #PERFORM OPTIMIZATION
    phi0 = .getPhi(thetaStart,nC,modTypes)
    tryCatch( { #This step may fail due to singularity
      suppressWarnings({
        foo = nlm(f=negloglikYphi, p=phi0,hessian=FALSE,iterlim=1000,steptol=steptol, progressbar=showProgressBar)#,print.level=2)
      })
      if(verbose && showProgressBar)  cat("\n") #ensure skipping line after computations
      nEvals = progcount + nEvals #update the number of evaluations (from progress)
      likVal = -foo$min #obtain estimated maximum likelihood value
      if(verbose) print(paste0("Maximum point found at: loglik=",round(likVal,2)))
  #    signif(.convBack(foo$est,nC, modTypes),2)

      #check if a better solution      
      if(likVal>maxL) {
        maxPhi = foo$estimate #store estimate
        maxL = likVal #store maximum
      }
    },error=function(e) print(e) )#,finally = {nITERinf <- nITERinf + 1} ) #end trycatch (update counter)
  } #end for each start values

############################OPTIM DONE#################
  validOptim = !is.na(maxPhi[1]) #check if valid optimization was found
  
  #Post-operations: Calculate 1) marker specific likelihood 2) Deconvolution 3) Model validation
  #1) Obtain loglik per marker:
  max_logliki = setNames(rep(NA,length(c$markerNames)),c$markerNames)
  DClist = list() #init as nothing (can also happen if full conditional)
  validTable = NULL
  maxHessian = NULL
  maxSigma <- matrix(NA,np,np)#Set as NA
  if(validOptim) {
    if(verbose) print("Optimizing done. Proceeding with post-calculations...")
    
    #Ensure that MxResult is correct order for the unknowns
    largePhi = 1000 #impute with this value if nan or inf (may affect hessian)
    maxPhi[is.nan(maxPhi)] = largePhi  #Handle that Mx-part of phi can have odd values on Mx-zero-boundaries
    thetahatFull = .convBack(maxPhi,nC, modTypes) #OBTAIN FULL PARAMETER LENGTH    
    thetahatUnknowns = .getThetaUnknowns(thetahatFull,nC,modTypes,inclLastMx=TRUE)
    thetahatUnknowns[1:nC] = .getMxValid(thetahatUnknowns[1:nC],nC,c$nU,c$hasKinship)
    maxPhi = .getPhi(thetahatUnknowns[-nC],nC,modTypes) #ensure a valid phi (convert back)

    #NEED TO CHECK IF ANY INFINITE ON PHI-estimate (FILL IN WITH LARGE NUMBER)
    if(nC>1) {   #Handle that Mx-part of phi can have odd values on Mx-zero-boundaries
      isMxRange = 1:(nC-1) #obtain range of Mx
      fillAsInf = is.infinite(maxPhi[isMxRange]) | is.nan(maxPhi[isMxRange])
      maxPhi[isMxRange[fillAsInf]] = largePhi #gives Mx=0
    }
    if(any(modTypes[-1])) { #check if any stutter model used:
      isStuttRange = nC + 1 + as.integer(modTypes[1]) + which(modTypes[-1]) #get Stutter range
      fillAsInf = is.infinite(maxPhi[isStuttRange]) | is.nan(maxPhi[isStuttRange])
      maxPhi[isStuttRange[fillAsInf]] = -largePhi #Gives StuttProp=0
    }
    #if(verbose && any(abs(validPhi-maxPhi) > 1e-6)) print("NOTE: Optimized Mx solution was reordred in order to be valid!")

    #CALCULATE COVARIANCE MATRIX:
    maxHessian = numDeriv::hessian(negloglikYphi,maxPhi,progressbar=FALSE)
    #maxHessian2 = Rdistance::secondDeriv(maxPhi,negloglikYphi,progressbar=FALSE)
    try({ #Avoid crash just in case maxHessian  is not valid
      maxSigma = MASS::ginv(maxHessian) #Use Pseudoinverse
    })
    thetahat2 = .convBack(maxPhi,nC, modTypes) #OBTAIN FULL PARAMETER LENGTH    obj$loglik(as.numeric(thetahat2) );   #repeat this to calc marker specific values for the MLE params
    obj$loglik(as.numeric(thetahat2) ); #repeat one more time to evaluate the "correct" likelihood
    #NOTE THAT HERE IT IS POSSIBLE TO CALL obj$calcGenoWeightsMax to get best precision.
    max_logliki = as.numeric( obj$logliki() ) #get marker specific likelihoods  (ensure it is a vector)
    names(max_logliki) = c$markerNames
    max_loglikiSUM = sum(max_logliki) #get sum 
    #negloglikYphi(phi=maxPhi,FALSE)  #equal
    if(max_loglikiSUM>maxL) maxL = max_loglikiSUM #set highest obtained

    #2) Obtain marginal probs of deconvolved profiles
    if(any(c$nUnknowns>0)) {
      DClist = obj$calcMargDC();   
      #sum(DClist[[1]]) # 8.711177e-13
      for(markerIdx in seq_along(DClist)) {
        DClist[[markerIdx]] = DClist[[markerIdx]]/rowSums(DClist[[markerIdx]]) #normalize
        #rownames(DClist[[m]]) = paste0("U",1:)
        genos = c$genoList[[markerIdx]] #obtain genotypes
        colnames(DClist[[markerIdx]]) = paste0(genos[,1],"/",genos[,2])
        if(c$nUnknowns[markerIdx]>0) rownames(DClist[[markerIdx]]) = paste0("U",seq_len(c$nUnknowns[markerIdx]))
      }
      names(DClist) = c$markerNames
    }
    #3) Obtain model validation
    validList = obj$calcValidMLE(); 
    for(markerIdx in seq_along(validList)) {
      UaPH = validList[[markerIdx]][1,]
      UaMAX = validList[[markerIdx]][2,]
      cumprobi = UaPH/UaMAX #obtaining cumulative probs for marker
      
      tableTemp = NULL #store in table
      for(aind in seq_len(c$nAlleles[markerIdx])) { #for each allele
        for(rind in seq_len(c$nRepMarkers[markerIdx])) { #for each replicate 
          cind = c$startIndMarker_nAllelesReps[markerIdx] + (aind-1)*c$nRepMarkers[markerIdx] + rind #get replicate-index
          alleleName  = c$alleleNames[aind + c$startIndMarker_nAlleles[markerIdx]] #obtain allele name
          newrow = data.frame( Marker=c$markerNames[markerIdx], Sample=c$repNames[ rind ], Allele=alleleName, Height=c$peaks[cind])
          tableTemp = rbind(tableTemp, newrow)
        }
      }
      tableTemp$ProbObs=as.numeric(cumprobi)
      tableTemp = tableTemp[tableTemp$Height>0,,drop=FALSE] #dont include dropouts
      validTable = rbind(validTable, tableTemp)
    }
  } #end if not valid optim
  
  #CLOSE after storing time
  time = ceiling(as.numeric(Sys.time() - start_time, units="secs")) #obtain time usage in seconds
  if(verbose) {
    print(paste0("Calculation time: ", .secondToTimeformat(time), " (HH:MM:SS)")) 
    print("-------------------------------------")
  } 
    
  obj$close() #free memory
   #.secondToTimeformat(time)

#########################################################
  #POST-PROCESSING... Calculate hessian
  fit = .getFittedParams(maxPhi,maxSigma,nC,modTypes)
  fit$maxHessian = maxHessian #store this as well
#  fit$thetaSE
#  fit$thetahat2
  
  detSIGMA = determinant(fit$thetaSigma)$mod[1] #calculate determinant
  if((is.na(detSIGMA) || is.infinite(detSIGMA)) && verbose) print("WARNING: The model is probably overstated, please reduce the model and fit it again!") 
  
  #laplace approx and handle if maxL is maximum size:
  if(abs(maxL)== .Machine$double.xmax) maxL = -Inf #its infinite
  logmargL <- 0.5*(np*log(2*pi)+detSIGMA) + maxL #get log-marginalized likelihood
  if(is.na(logmargL)) logmargL = -Inf #special handle here
  
  nUnknown = c$nU - as.integer(c$hasKinship) #number of Mx params to be restricted
  if(nUnknown>1) { #adjust if more than 1 unknown 
    logmargL <- lfactorial(nUnknown) + logmargL #log(factorial(nU)) + logmargL #get adjusting for symmetry of unknowns
  }
  #prepare return variable
  fit <- c(fit,list(loglik=maxL,logmargL=logmargL, logliki=max_logliki, nEvals=nEvals, time=time)) #append
  model <- list(nC=nC,nU=c$nU,samples=samples,popFreq=popFreq,refData=refData,condOrder=condOrder,knownRef=knownRef,BWS=BWS,FWS=FWS,DEG=DEG,prC=pC, AT=AT,fst=fst,
                lambda=lambda,kit=kit,minF=minF,normalize=normalize,knownRel=knownRel,ibd=ibd,priorBWS=priorBWS,priorFWS=priorFWS, adjQbp=adjQbp)
  ret <- list(fit=fit,model=model,nDone=nDone,delta=delta,steptol=steptol,seed=seed,prepareC=c, difftol=difftol,resttol=resttol,thetaPresearch=thetaPresearch,restprop=restprop) 
  ret$model$threshT = AT #make a copy (used in plotTop functions)
  ret$modelType = modTypes #obtain model type directly)
  ret$DCmargList = DClist #attach the DC list 
  ret$MLEvalidTable = validTable #attach the MLE validation list 
  return(ret)
}