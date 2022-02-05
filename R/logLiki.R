#' @title logLiki
#' @author Oyvind Bleka
#' @description logLiki returns the likelihood for each marker based on the MLE fit
#'
#' @param mlefit Fitted object using contLikMLE
#' @param verbose Whether printing deconvolution progress. Default is TRUE.
#' @param maxThreads Maximum number of threads to be executed by the parallelization
#' @return ret A vector with log-likelihood-values for each locus for given model
#' @export

logLiki <- function(mlefit,verbose=FALSE,maxThreads=32){
  
  thhat <- mlefit$fit$thetahat #condition on mle parameter
  model <- mlefit$model #take out assumed model with given data
  locs <- names(model$popFreq) #get loci names to evaluate
  c <- mlefit$prepareC #returned from prepareC
  nC = c$nC #number of contributors
  usedeg <- !is.null(model$kit) #check for degradation
  
  #Prepare fixed params (create vector if not already):
  nM = c$nM #number of markers to evaluate
  if(length(locs)!=nM) stop("Number of loci from prepareC and popFreq in mlefit was not the same!")
  if(any(is.na(thhat))) return(rep(NA,nM)) #return empty values (NA will give NA per markers)
  
  ATv = model$threshT
  pCv = model$prC
  lambdav = model$lambda
  fstv = model$fst

  #Prepare estimated params:
  mixprop = head(thhat,nC-1)# fitted mix proportions
  mu1 = thhat[nC] #fitted PHexp
  sig1 = thhat[nC+1] #fitted PHvar
  
  #Obtain remaining params
  beta1 = 1 #set default values
  xiB = model$xi #set default values
  xiF = model$xiFW #set default values
  indsel = nC + 2  #get index to insert estimated variable form param vector(thhat)
  if(!is.null(model$kit)) {
    beta1 = thhat[indsel] #get degrad param
    indsel = indsel + 1 #update index
  }
  if(is.null(xiB)) {
    xiB = thhat[indsel] #get BWstutter param
    indsel = indsel + 1 #update index
  }
  if(is.null(xiF)) xiF = thhat[indsel] #get FWstutter param
  
  logLikv = rep(NA,nM) #obtain loglik per marker
  names(logLikv) = locs #name marker names
  
  #Step 1) Calculate L(E|thetahat) for each marker
  startIndPS = 0; #start marker index for potential stutters (own vectors)
  startIndMarker1=0;#start marker index (1 rep)
  startIndMarker2=0;#start marker index (nRep[m] reps)

  progcount  = 1 #progress counter
  if(verbose) progbar <- txtProgressBar(min = 0, max =nM, style = 3) #create progress bar
  for(m in 1:nM) { #extract info in c relevant for each markers:
    ind1 = startIndMarker1 + 1:c$nA[m] #get index of 1 repitition
    ind2 = startIndMarker2 + 1:(c$nA[m]*c$nRep[m]) #get index of 1 repitition
    indPS = startIndPS +  1:c$nPS[m] #get index of potential stutters
    if(c$nPS[m]==0) indPS = numeric()
	
  	genoinds = (nC*(m-1)+1):(nC*m) #get index of known genotype index (knownGind and relGind)    
    logLikv[m] = .C("loglikgammaC",as.numeric(0),as.integer(nC),as.integer(c$NOK[m]),as.integer(c$knownGind[genoinds] ),as.numeric(mixprop),as.numeric(mu1),as.numeric(sig1),as.numeric(beta1),as.numeric(xiB),as.numeric(xiF),as.numeric(ATv[m]),as.numeric(pCv[m]),as.numeric(lambdav[m]),as.numeric(fstv[m]),c$nReps[m],as.integer(1),c$nA[m],c$YvecLong[ind2],c$FvecLong[ind1],c$nTypedLong[m],c$maTypedLong[ind1],c$basepairLong[ind1],c$BWvecLong[ind1],c$FWvecLong[ind1],c$nPS[m],c$BWPvecLong[indPS],c$FWPvecLong[indPS],as.integer(maxThreads),as.integer(0),as.integer(c$anyRel),as.integer(c$relGind[genoinds]),c$ibdLong,PACKAGE="euroformix")[[1]]

    #Update indices for next marker:
    startIndMarker1 = startIndMarker1 + c$nA[m] #get start position of marker m+1 (1 rep)
    startIndMarker2 = startIndMarker2 + c$nA[m]*c$nRep[m]; #get start position of marker m+1 (nRep), used only for Peaks only
    startIndPS = startIndPS + c$nPS[m]; #add number of potential stutters
    
    if(verbose) { #only show progressbar if verbose (update)
      setTxtProgressBar(progbar,progcount)
      progcount <- progcount + 1
    } 
    
  }    
  
 return( logLikv ) #simply returns  
}



