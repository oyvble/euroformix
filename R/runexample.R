#' @title runexample
#' @author Oyvind Bleka
#' @description Calculate LR based on the inbuilt methods in euroformix for a specific setting
#' @details Using following inbuilt methods:  sample_tableToList,prepareData,contLikMLE,logLiki,validMLEmodel,deconvolve,contLikMCMC,contLikINT, calcGjoint,plotEPG2
#' @param NOC Number of contributors
#' @param popfn Filename of population frequency file (must contain EFM compatible text file)
#' @param evidfn Filename of evidence(s) profile file (must contain EFM compatible text file)
#' @param reffn Filename of reference(s) profile file (must contain EFM compatible text file)
#' @param POIind Index of the references being Person of interest (POI), only Hp
#' @param condInd Index of conditional references (both Hp and Hd)
#' @param kit Kitname (shortname of kit obtained by getKit()). Used to model degradation and visualize EPG. 
#' @param modelDegrad Boolean of whether degradation should be modeled (requires kit to be specified)
#' @param modelBWstutt Boolean of whether backward stutter should be modeled
#' @param modelFWstutt Boolean of whether forward stutter should be modeled
#' @param findOptimalModel Boolean of whether to find the optimal model (using AIC criterion)
#' @param fst The co-ancestry coefficient (theta-correction). Can be a vector (must be number of marker long)
#' @param lambda Parameter in modeled peak height shifted exponential model. Can be a vector (must be number of marker long)
#' @param prC Allele drop-in probability. Can be a vector (must be number of marker long)
#' @param threshT The analytical/detection threshold given. Can be a vector (must be number of marker long, and named with marker names to plot)
#' @param alpha The significance level used for the envelope test (validMLEmodel). Default is 0.01
#' @param minF Minimum frequency used in the analysis (can be specified, otherwise minimum observed is used by default)
#' @param seed Seed used in contLikMLE/contLikMCMC)
#' @param nDone Number of required optimizations in contLikMLE 
#' @param mcmcIter Number of iterations used in contLikMCMC
#' @param mcmcDelta Scaling of variance used for proposal in contLikMCMC (calibrated based on mcmc$accrat object)
#' @param intMaxEval max number of evalutions in integration (contLikINT)
#' @param intRelTol relative tolerance for integration  (contLikINT)
#' @param verbose Boolean whether printing data (EPG) and progress
#' @param maxThreads Maximum number of threads to be executed by the parallelization
#' @return A list with calculated result elements
#' @examples
#' \dontrun{
#' kit = "ESX17"
#' popfn = paste(path.package("euroformix"),"FreqDatabases",
#'  paste0(kit,"_Norway.csv"),sep=.Platform$file.sep)
#' evidfn = paste(path.package("euroformix"),"examples",paste0(kit,"_3p.csv"),
#'  sep=.Platform$file.sep)
#' reffn = paste(path.package("euroformix"),"examples",paste0(kit,"_refs.csv"),
#'  sep=.Platform$file.sep)
#' resList = runexample(NOC=3,popfn,evidfn,reffn,POIind=1,condInd=2,kit=kit,
#'  modelDegrad=TRUE,modelBWstutt=TRUE,verbose=TRUE)
#' }
#' @export
runexample = function(NOC,popfn,evidfn,reffn,POIind=1,condInd=NULL,kit=NULL, modelDegrad=FALSE,modelBWstutt=FALSE,modelFWstutt=FALSE,findOptimalModel=TRUE,fst=0.01,lambda=0.01,prC=0.05,threshT=50,alpha=0.01,minF=NULL,seed=1,nDone=2,mcmcIter=500,mcmcDelta=2,intRelTol=0.1,intMaxEval=1000,verbose=TRUE,maxThreads=32) { 
  LUSsymbol = "_" #LUS symbold
  
  if(modelDegrad && is.null(kit)) stop("Please specify a valid kit to model degradation")  
  
  #Specify optional model setup:
  kit0 = kit #store kit input
  if(!modelDegrad) kit = NULL #degradation model turned off when kit=NULL
  xi <- xiFW <- 0 #default is no stutters
  if(modelBWstutt) xi = NULL
  if(modelFWstutt) xiFW = NULL
  
  #Obtain data:
  popFreq = freqImport(popfn)[[1]] #obtain list with population frequencies
  samples = sample_tableToList(tableReader(evidfn))#,threshT)
  refs = sample_tableToList(tableReader(reffn))
  numRefs = length(refs) #get number of refs
  if( numRefs==0 ) stop("No references provided. Please specify a valid reference file.")

  #Determine data type:
  isMPS = FALSE
  isLUS = all(unlist(sapply(samples,function(x)  sapply(x,function(y) all(grepl(LUSsymbol,y$adata),na.rm=TRUE)) )))  
  if(!isLUS) isMPS = all(unlist(sapply(samples,function(x)  sapply(x,function(y) all( is.numeric(y$adata) )))))
  if(is.null(kit0)) isMPS = TRUE    #print("Please specify kit in order to show EPG. Otherwise default MPS format is shown!")
  isEPG = !isMPS && !isLUS  #is EPG if not LUS OR MPS
  
  
  #Set up hypothesis (contributors)
  cond = rep(0,numRefs) #init condition vector (used to decide which to condition on)
  if(!is.null(condInd)) cond[condInd] = 1:length(condInd) #set contributors
  condhp <- condhd <- cond #copy variable
  condhp[POIind] = max(condhd) + 1 #add POI contributor
  knownRefhd = POIind #known reference under Hd

  if( sum(condhp>0) > NOC ) stop("The number of contributors were less than the number of conditionals. Please increase the number of contributors.")
    
  #Show data:
  if(verbose) {
    if(isEPG) {
      plotEPG2(samples,kit0,refs, AT=threshT) 
    } else {
      plotMPS2(samples,refs, AT=threshT) 
    }
  }
  dat = prepareData(samples,refs,popFreq,threshT=threshT,minF=minF, normalize = FALSE) #obtain data to use for analysis
 
  #obtain random match probability for each markeres
  rmp = calcRMPfst(dat,POIind,condInd,fst=fst ) #non-scaled values
  LRupper = -sum(log10(rmp)) #maximum attainable LR (log10 scale)

  modelSelTable = NULL
  if(findOptimalModel) {
    
    #Running through all configurations (uses contLikMLE/validMLEmodel functions):
    modelDegrad = unique(c(FALSE,modelDegrad))
    modelBWstutt = unique(c(FALSE,modelBWstutt))
    modelFWstutt = unique(c(FALSE,modelFWstutt))
    minNOC = sum(condhp>0)
    searchList = contLikSearch( NOC= minNOC:NOC, modelDegrad,modelBWstutt,modelFWstutt,samples=samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=condhd,knownRefPOI=POIind, prC=prC,threshT=threshT,fst=fst,lambda=lambda,kit=kit,nDone=nDone,seed=seed,verbose=verbose,maxThreads=maxThreads,alpha=alpha)
    adjAIC = searchList$outtable[,2] #obtain crietion
    optimInd = which.max(adjAIC)[1] #get index of optimal model. Use simpler model if "Tie"
      #set to optimal model:
    mlehp = searchList$hpfitList[[optimInd]]
    mlehd = searchList$hdfitList[[optimInd]]
    modelSelTable = searchList$outtable #get model selection table
    
  } else { #otherwise use current model
    #if(verbose) print("Performing Maximum likelihood estimation (optimize)...")
    mlehp = contLikMLE(nC=NOC,samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=condhp,xi=xi,prC=prC,nDone=nDone,threshT=threshT,fst=fst,lambda=lambda,kit=kit,xiFW=xiFW,seed=seed,verbose=verbose,maxThreads=maxThreads)
    mlehd = contLikMLE(nC=NOC,samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=condhd,xi=xi,prC=prC,nDone=nDone,threshT=threshT,fst=fst,lambda=lambda,kit=kit,xiFW=xiFW,seed=seed,verbose=verbose,maxThreads=maxThreads,knownRef = knownRefhd)
    }
  #MLE VALIDATION:
  if(verbose) print("Performing MLE validations...")
  validhp = validMLEmodel(mlehp,kit0,plottitle = "Hp",maxThreads=maxThreads,alpha=alpha)
  validhd = validMLEmodel(mlehd,kit0,plottitle = "Hd",maxThreads=maxThreads,alpha=alpha)
  
  #calculate logLik per markers:
  if(verbose) print("Calculating logLik per markers...")
  hpv = logLiki(mlefit=mlehp,verbose=verbose,maxThreads=maxThreads)
  hdv = logLiki(mlefit=mlehd,verbose=verbose,maxThreads=maxThreads)
  all.equal(sum(hpv),mlehp$fit$loglik) #check that its same as joint loglik
  all.equal(sum(hdv),mlehd$fit$loglik) #check that its same as joint loglik
  LRmle = (mlehp$fit$loglik - mlehd$fit$loglik)/log(10) #obtain LR based on mle (log10 scale)
  LRmarker = (hpv - hdv)/log(10) #obtain LR per marker  (log10 scale)
      
  #Perform DC:
  if(verbose) print("Performing Deconvolution...")
  DChp = deconvolve(mlefit=mlehp,verbose=verbose)
  DChd = deconvolve(mlefit=mlehd,verbose=verbose)
#  DChd$table2
  
  #Show expected PHs with data:
  if(verbose) {
    if(isEPG) {
      plotTopEPG2(mlehp,DChp,kit0) #show fitted PH under Hp
      plotTopEPG2(mlehd,DChd,kit0) #show fitted PH under Hd
    } else {
      plotTopMPS2(mlehp,DChp) #MLEobj=mlehp;DCobj=DChp;grpsymbol="_";locYmax=TRUE;options=NULL
      plotTopMPS2(mlehd,DChd) #show fitted PH under Hd
    }
  }
  
  #Perform MCMC:
  if(verbose) print("Performing MCMC simulations...")
  mcmchp = contLikMCMC(mlehp,niter=mcmcIter,delta=mcmcDelta,seed=seed,verbose=verbose,maxThreads=maxThreads)
  if(verbose) print(paste0("Acceptance rate under Hp=",mcmchp$accrat," (should be around 0.2)"))
  mcmchd = contLikMCMC(mlehd,niter=mcmcIter,delta=mcmcDelta,seed=seed + 999,verbose=verbose,maxThreads=maxThreads)
  if(verbose) print(paste0("Acceptance rate under Hd=",mcmchd$accrat," (should be around 0.2)"))
  
  validMCMC(mcmchp) #diagnostic of mcmc 
  validMCMC(mcmchd) #diagnostic of mccm
  LRdistr = (mcmchp$postlogL - mcmchd$postlogL)/log(10) #get LR distirbution
  LRmarg = log10(mcmchp$margL/mcmchd$margL) #estimated marginalised
  #plot(density(LRdistr));abline(v=LRmle,lty=2)
  LRcons = quantile(LRdistr,0.05) #extract 5% quantile  (conservative LR)
  
  #Marginalised likelihood estimation (optimize)
  LRlaplace = (mlehp$fit$logmargL - mlehd$fit$logmargL)/log(10) #obtain LR estimate based on Laplace approx

  if(verbose) print("Performing Numerical integrations...")
  #Under Hp:
  rng = apply(mcmchp$posttheta,2,range) #get range of parameters (under Hp)
  scale= -mlehp$fit$loglik  #scale to adjust loglikelihood
  intHp = contLikINT(lower=rng[1,],upper=rng[2,], nC=NOC,samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=condhp,xi=xi,prC=prC,threshT=threshT,fst=fst,lambda=lambda,kit=kit,reltol=intRelTol,scale=scale,xiFW=xiFW,maxEval=intMaxEval,verbose=verbose,maxThreads=maxThreads)
  loglikINThp = log(intHp$margL)-scale #obtain integral

  #Under Hd:
  rng = apply(mcmchd$posttheta,2,range) #get range of parameters (under Hp)
  scale= -mlehd$fit$loglik  #scale to adjust loglikelihood
  intHd = contLikINT(lower=rng[1,],upper=rng[2,], nC=NOC,samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=condhd,xi=xi,prC=prC,threshT=threshT,fst=fst,lambda=lambda,kit=kit,reltol=intRelTol,scale=scale,xiFW=xiFW,maxEval=intMaxEval,verbose=verbose,maxThreads=maxThreads,knownRef = knownRefhd)
  loglikINThd = log(intHd$margL)-scale #obtain integral
  
  LRint = (loglikINThp - loglikINThd)/log(10)
  
  return( list(LRupper=LRupper,LRmle=LRmle,LRcons=LRcons,LRint=LRint,LRlaplace=LRlaplace,LRmarg=LRmarg,DChp=DChp$table2,DChd=DChd$table2,validhp=validhp,validhd=validhd,LRmarker = LRmarker,modelSelTable=modelSelTable) )
    
}
