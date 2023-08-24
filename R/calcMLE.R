
#' @title calcMLE
#' @author Oyvind Bleka
#' @description Optimizes the likelihood function of the DNA mixture model 
#' @param nC Number of contributors in model.
#' @param samples A List with samples which for each samples has locus-list elements with list elements adata and hdata. 'adata' is a qualitative (allele) data vector and 'hdata' is a quantitative (peak heights) data vector.
#' @param popFreq A list of allele frequencies for a given population.
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
#' @param steptol Argument used in the nlm function for faster return from the optimization (tradeoff is lower accuracy).
#' @param nDone Number of optimizations required providing equivalent results (same logLik value obtained)
#' @param delta Scaling of variation of normal distribution when drawing random startpoints. Default is 1.
#' @param difftol Tolerance for being exact in log-likelihood value (relevant when nDone>1)
#' @param seed The user can set seed if wanted 
#' @param verbose Whether printing limits to integrate over. Printing progress if maxEval>0. Default is TRUE.
#' @param priorBWS Prior function for BWS-parameter. Flat prior on [0,1] is default.
#' @param priorFWS Prior function for FWS-parameter. Flat prior on [0,1] is default.
#' @param maxThreads Maximum number of threads to be executed by the parallelization
#' @param adjQbp Indicate whether fragmenth length of Q-allele is based on averaged weighted with frequencies
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
#' mlefit = calcMLE(3,samples,popFreq,refData,condOrder,kit=kit)
#' plotTopEPG2(mlefit)
#' }

calcMLE = function(nC,samples,popFreq, refData=NULL, condOrder = NULL, knownRef = NULL, kit=NULL,DEG=TRUE,BWS=TRUE,FWS=TRUE,
  AT=50,pC=0.05,lambda=0.01,fst=0,knownRel=NULL,ibd=NULL,minF=NULL,normalize=TRUE,steptol = 1e-4,nDone=3,delta=1,difftol=0.01,
  seed=NULL,verbose=FALSE,priorBWS=NULL,priorFWS=NULL, maxThreads=0, adjQbp = FALSE) {

  start_time <- Sys.time()
  if(!is.null(seed)) set.seed(seed) #set seed if provided

  #Preparing variables to fill into C++ structure
  c = prepareC(nC,samples,popFreq, refData, condOrder, knownRef, kit,BWS,FWS,AT,pC,lambda,fst,knownRel,ibd,minF,normalize, adjQbp)
  if(DEG && !c$useDEG) {
    print("Degradation had to be turned off since kitinfo was not found for selected kit.")
    DEG = FALSE
  }
  modTypes = setNames(c(DEG,BWS,FWS),c("DEG","BWS","FWS")) #obtain model types
  #PREPARE FEEDING DATA
  if(verbose) print("Carrying out preparations for optimizing...")
  mod = Rcpp::Module( "mod",PACKAGE="euroformix" ) #load module
  obj = methods::new(mod$ExposedClass) #create object of class
  
  #Step 1: insert data to exposed class (filldata)
  obj$filldata(c$nStutterModels,c$nMarkers,c$nRepMarkers,c$nAlleles,c$startIndMarker_nAlleles,c$startIndMarker_nAllelesReps,c$peaks,c$freqs,c$dropinWeight, c$nTyped, c$maTyped, c$basepair,
               c$BWfrom, c$FWfrom, c$BWto, c$FWto, c$nPotStutters, c$startIndMarker_nAllelesTot, c$QalleleIndex, c$dropinProb, c$fst, c$AT, c$NOK, c$knownGind, c$relGind, c$ibd, as.integer(maxThreads)) 
  
  #Step 2: Indexing large matrix (doIndex)
  prepTime=system.time({
    obj$prepare(as.integer(nC))
  })[3]
  #sum((nA+nPS)*(nA*(nA+1)/2)^nC) #number of evaluations
  if(verbose) print(paste0("Prep. done and took ",round(prepTime),"s. Start optimizing..."))
  
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
  ncond = sum(condOrder>0) #numbers to conidtion on
  nU = nC-ncond
  np = (nC+1) + DEG + BWS + FWS  #obtain number of parameters
  
  if(!DEG) meanbp = NULL
  #obtain expected startpoints of (PHexp,PHvar,DegSlope)
  th0 <- fitgammamodel(y=sumY,x=meanbp,niter=10,delta=delta,offset=0,scale=1) 

  #function for calling on C-function: Must convert "real domain" values back to model params 
  negloglikYphi <- function(phi,progressbar=TRUE) { #assumed order: mixprop(1:C-1),mu,sigma,beta,xi
    param = .convBack(phi,nC, modTypes)
    loglik = obj$loglik(as.numeric(param) ) 
    paramStutters = param[length(param) + c(-1,0)] #obtain stutter params (last two)
    if(!is.null(priorBWS) && paramStutters[1]>0) loglik = loglik + log( priorBWS(paramStutters[1]) )
    if(!is.null(priorFWS) && paramStutters[2]>0) loglik = loglik + log( priorFWS(paramStutters[2]) )
    
    progcount <<- progcount + 1 #update counter
    if(verbose && progressbar) setTxtProgressBar(progbar,progcount)  #only show progressbar if verbose(and if decided to show)
    return(-loglik) #weight with prior of stutter.
  }  

  #Strategy to obtain global maximum: Assure that maximum logLik (maxL) is obtained in 'nDone' optimizations 
  maxITERinf <- 30 #number of possible times to be INF or not valid optimum before any acceptance
  nITERinf <- 0 #number of times beeing INF (invalid value)
  maxL <- -Inf #value of maximum obtained loglik
  nOK <- 0 #number of times for reaching largest previously seen optimum
  maxIterProgress = (np-2)*100 #maximum iterations in progressbar
  nEvals = 0 #number of total evaluations (calls to likelihood function)
  
  suppressWarnings({
    while(nOK<nDone) {
      #Obtain random start values for parameters
      p0 <- .paramrandomizer(th0,nC,modTypes,delta,ncond,c$hasKinship)#,T) #generate random start value on Real (don't need to)
      #prettyNum(convBack(p0,nC, modTypes))
      
      progcount  = 1 #progress counter for optimization (reset for each sucessful optimization)
      timeOneCall = system.time({ #estimate the time for calling the likelihood function one time 
        likval <- negloglikYphi(phi=p0,FALSE)   #check if start value was accepted
      })[3] #obtain time in seconds
      
      if( is.infinite(likval) ) { #if it was infinite (invalid)
        nITERinf = nITERinf + 1	 
        
      } else {
        if(verbose) {
          showProgressBar = FALSE #show progress bar only if calculatations take more than 10 seconds
          expectedTimeProgress0 = timeOneCall*maxIterProgress
          if(expectedTimeProgress0 > 10) showProgressBar = TRUE #show progress if upper time >10s
          if(showProgressBar) {
            expectedTimeProgress = .secondToTimeformat(expectedTimeProgress0)
            cat(paste0("\nExpected (upper) time is ", expectedTimeProgress, " (HH:MM:SS):\n")) 
            #\n---------------------------------------------------
          }
        }
        
        #if(verbose) progbar <- tcltk::tkProgressBar(min = 0, max = iterscale*maxIter,width = 300) #create progress bar
        if(verbose && showProgressBar) progbar <- txtProgressBar(min = 0, max = maxIterProgress, style = 3) #create progress bar
        
        tryCatch( {
          foo <- nlm(f=negloglikYphi, p=p0,hessian=TRUE,iterlim=500,steptol=steptol, progressbar=showProgressBar)#,print.level=2)
          nEvals = progcount + nEvals #update the number of evaluations (from progress)
          Sigma <- solve(foo$hessian)
          
          if(all(diag(Sigma)>=0) && foo$iterations>2) { #} && foo$code%in%c(1,2)) { #REQUIREMENT FOR BEING ACCEPTED
            nITERinf <- 0 #reset INF if accepted
            likval <- -foo$min #obtain local maximum
            
            #was the maximum (approx) equal the prev: Using decimal numbers as difference 
            isEqual = !is.infinite(maxL) && abs(likval-maxL)< difftol  #all.equal(likval,maxL, tolerance = 1e-2) # # #was the maximum (approx) equal the prev?
            
            if(isEqual) { 
              if(verbose) {
                if(showProgressBar) cat("\n") #skip line
                print(paste0("Equal maximum found: loglik=",likval))
              }
              nOK = nOK + 1 #add counter by 1
            } else {  #if values were different
              if(likval>maxL) { #if new value is better
                nOK = 1 #first accepted optimization found
                maxL <- likval #maximized likelihood
                maxPhi <- foo$est #set as topfoo     
                maxSigma <- Sigma 
                if(verbose) {
                  if(showProgressBar)  cat("\n") #skip line
                  print(paste0("New maximum at loglik=",likval))
                } 
                
                # if(verbose) cat(paste0("maxPhi=",paste0(maxPhi,collapse=","))) #removed in v2.0
              } else {
                if(verbose) {
                  if(showProgressBar)  cat("\n") #skip line
                  print(paste0("Local (non-global) maximum found at logLik=",likval))
                }
              }
            } 
            if(verbose) print(paste0(" (",nOK,"/",nDone,") optimizations done"))
          } else { #NOT ACCEPTED
            nITERinf <- nITERinf + 1 
          }
        },error=function(e) e,finally = {nITERinf <- nITERinf + 1} ) #end trycatch (update counter)
        
        if(verbose && showProgressBar)  cat("\n") #ensure skipping line after computations
      } #end if else 
      
      if(nOK==0 && nITERinf>maxITERinf) {
        nOK <- nDone #finish loop
        maxL <- -Inf #maximized likelihood
        maxPhi <- rep(NA,np) #Set as NA
        maxSigma <- matrix(NA,np,np)#Set as NA
        break #stop loop if too many iterations  
      }
    } #end while loop
  })    
  #Obtain params:
  thetahat2 = .convBack(maxPhi,nC, modTypes) #OBTAIN FULL PARAMETER LENGTH

  #Obtain loglik per marker:
  obj$loglik(as.numeric(thetahat2) );   #repeat one more time to calc log-specific
  max_logliki = as.numeric( obj$logliki() ) #get marker specific likelihoods  (ensure it is a vector)
  names(max_logliki) = c$markerNames
  obj$close() #free memory
   #.secondToTimeformat(time)

  #POST-PROCESSING (COVARIANCE MATRIX)
  if(any(!modTypes)) thetahat2 = thetahat2[-(nC+2+which(!modTypes))] #remove non-modelled values
  thetahat = thetahat2[-nC] #also remove MxProp (nC)

  #OBTAIN COVARIANCE MATRIX OF PARAMS:
  J <- .calcJacobian(maxPhi,thetahat,nC) #obtain jacobian matrix:
  Sigma <- (t(J)%*%maxSigma%*%J) #this is correct covariance of thetahat. Observed hessian is used
  
  #get extended Sigma (all parameters)
  Sigma2 <- matrix(NA,nrow=np+1,ncol=np+1) #extended covariance matrix also including mx[nC]
  Sigma2[nC:np+1,nC:np+1] <- Sigma[nC:np,nC:np] 
  Sigma2[nC:np+1,nC:np+1] <- Sigma[nC:np,nC:np] 
  if(nC>1) {
    Sigma2[nC:np+1,1:(nC-1)] <- Sigma[nC:np,1:(nC-1)] 
    Sigma2[1:(nC-1),1:(nC-1)] <- Sigma[1:(nC-1),1:(nC-1)] 
    Sigma2[1:(nC-1),nC:np+1] <- Sigma[1:(nC-1),nC:np] 
    Sigma2[nC,nC] <- sum(Sigma[1:(nC-1),1:(nC-1)])
    for(k in (1:(np+1))[-nC]) {
      Sigma2[nC,k] <- Sigma2[k,nC] <- -sum(Sigma[1:(nC-1),k-sum(k>nC)]) 
    }
  } else {
    Sigma2[1,1:(np+1)] <- Sigma2[1:(np+1),1] <- 0 #no uncertainty
  }
  
  #Standard error for theta:
  thetaSE <- sqrt(diag(Sigma2))
  
  mxName <- "Mix-prop. C" #"mx"  
  thetanames0 <- c("P.H.expectation","P.H.variability")
  phinames <- paste0("log(",thetanames0,")")
  if(nC>1) {
    phinames  <- c(paste0("nu",1:(nC-1)),phinames)
    thetanames <- c(paste0(mxName,1:(nC-1)),thetanames0)
  } else {
    thetanames <- thetanames0
  }
  thetanames2 <- c(paste0(mxName,1:nC),thetanames0)
  
  othernames <- character()
  if(DEG) othernames = c(othernames,"Degrad. slope")
  if(BWS) othernames = c(othernames,"BWstutt-prop.")  
  if(FWS) othernames = c(othernames,"FWstutt-prop.")
  if(length(othernames)>0) phinames <- c(phinames,paste0("logit(",othernames,")") )
  thetanames <- c(thetanames, othernames )
  thetanames2 <- c(thetanames2, othernames )
  
  colnames(maxSigma) <- rownames(maxSigma) <- phinames 
  colnames(Sigma) <- rownames(Sigma) <- thetanames
  colnames(Sigma2) <- rownames(Sigma2) <- thetanames2
  names(maxPhi) <- phinames
  names(thetahat) <- thetanames
  names(thetahat2) <- thetanames2
  names(thetaSE) <- thetanames2
  
  #laplace approx:
  logmargL <- 0.5*(np*log(2*pi)+determinant(Sigma)$mod[1]) + maxL #get log-marginalized likelihood
  #nU <- nC #number of contributors
  if(nU>1) { #if more than 1 unknown 
    logmargL <- lgamma(nU+1) + logmargL #log(factorial(nU)) + logmargL #get adjusting for symmetry of unknowns
  }
  time = ceiling(as.numeric(Sys.time() - start_time, units="secs")) #obtain time usage in seconds
  
  fit <- list(phihat=maxPhi,thetahat=thetahat,thetahat2=thetahat2,phiSigma=maxSigma,thetaSigma=Sigma,thetaSigma2=Sigma2,loglik=maxL,thetaSE=thetaSE,logmargL=logmargL, logliki=max_logliki, nEvals=nEvals, time=time)
  model <- list(nC=nC,nU=nU,samples=samples,popFreq=popFreq,refData=refData,condOrder=condOrder,knownRef=knownRef,BWS=BWS,FWS=FWS,DEG=DEG,prC=pC, AT=AT,fst=fst,
                lambda=lambda,kit=kit,minF=minF,normalize=normalize,knownRel=knownRel,ibd=ibd,priorBWS=priorBWS,priorFWS=priorFWS, adjQbp=adjQbp)
  ret <- list(fit=fit,model=model,nDone=nDone,delta=delta,steptol=steptol,seed=seed,prepareC=c, difftol=difftol) 
  ret$model$threshT = AT #make a copy (used in plotTop functions)
  ret$modelType = modTypes #obtain model type directly)
  ret$maxThreads = as.integer(maxThreads)
  return(ret)
}