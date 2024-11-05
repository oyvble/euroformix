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
#' @param modelDegrad Whether degradation should be modeled (requires kit to be specified)
#' @param modelBWstutt Whether backward stutter should be modeled
#' @param modelFWstutt Whether forward stutter should be modeled
#' @param findOptimalModel Whether to find the optimal model (using AIC criterion)
#' @param fst The co-ancestry coefficient (theta-correction). Can be a vector (must be number of marker long)
#' @param lambda Parameter in modeled peak height shifted exponential model. Can be a vector (must be number of marker long)
#' @param prC Allele drop-in probability. Can be a vector (must be number of marker long)
#' @param threshT The analytical/detection threshold given. Can be a vector (must be number of marker long, and named with marker names to plot)
#' @param alpha The significance level used for the envelope test (validMLEmodel). Default is 0.01
#' @param seed Seed used in contLikMLE/contLikMCMC)
#' @param nDone Number of required optimizations in contLikMLE 
#' @param mcmcIter Number of iterations used in contLikMCMC
#' @param nTippets Number of iterations used for non-contributor/tippet analysis 
#' Scaling of variance used for proposal in contLikMCMC (calibrated based on mcmc$accrat object)
#' @param verbose Boolean whether printing data (EPG) and progress
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
#library(euroformix);NOC=2;POIind=1;condInd=2;modelDegrad=TRUE;modelBWstutt=TRUE;modelFWstutt=FALSE;verbose=FALSE
#kit = "ESX17";popfn = paste(path.package("euroformix"),"FreqDatabases",paste0(kit,"_Norway.csv"),sep=.Platform$file.sep);evidfn = paste(path.package("euroformix"),"examples",paste0(kit,"_3p.csv"),sep=.Platform$file.sep);reffn = paste(path.package("euroformix"),"examples",paste0(kit,"_refs.csv"),sep=.Platform$file.sep);
#fst=0.01;lambda=0.01;prC=0.05;threshT=50;alpha=0.01;minF=NULL;normalize=TRUE;seed=1;nDone=2;mcmcIter=500;mcmcDelta=2;intRelTol=0.1;intMaxEval=1000;findOptimalModel=FALSE
runexample = function(NOC,popfn,evidfn,reffn,POIind=1,condInd=NULL,kit=NULL, modelDegrad=FALSE,modelBWstutt=FALSE,modelFWstutt=FALSE,findOptimalModel=TRUE,fst=0.01,lambda=0.01,prC=0.05,threshT=50,alpha=0.01,seed=1,nDone=2,mcmcIter=1000,nTippets=10,verbose=TRUE) { 
  LUSsymbol = "_" #LUS symbold
  if(modelDegrad && is.null(kit)) stop("Please specify a valid kit to model degradation")  
  
  #Obtain data:
  popFreq = freqImport(popfn)[[1]] #obtain list with population frequencies
  samples = sample_tableToList(tableReader(evidfn))#,threshT)
  refData = sample_tableToList(tableReader(reffn))
  numRefs = length(refData) #get number of refData
  if( numRefs==0 ) stop("No references provided. Please specify a valid reference file.")

  #Determine data type:
  isMPS = FALSE
  isLUS = all(unlist(sapply(samples,function(x)  sapply(x,function(y) all(grepl(LUSsymbol,y$adata),na.rm=TRUE)) )))  
  if(!isLUS) isMPS = all(unlist(sapply(samples,function(x)  sapply(x,function(y) all( is.numeric(y$adata) )))))
  if(is.null(kit)) isMPS = TRUE    #print("Please specify kit in order to show EPG. Otherwise default MPS format is shown!")
  isEPG = !isMPS && !isLUS  #is EPG if not LUS OR MPS
  
  #Set up hypothesis (contributors)
  cond = rep(0,numRefs) #init condition vector (used to decide which to condition on)
  if(!is.null(condInd)) cond[condInd] = 1:length(condInd) #set contributors
  condhp <- condhd <- cond #copy variable
  condhp[POIind] = max(condhd) + 1 #add POI contributor
  knownRefhd = POIind #known reference under Hd

  if( sum(condhp>0) > NOC ) stop("The number of contributors were less than the number of conditionals. Please increase the number of contributors.")
#  mlefitHp = calcMLE(NOC,samples,popFreq,refData,condhp,NULL,       kit,modelDegrad,modelBWstutt,modelFWstutt,threshT,prC,lambda,fst,nDone=nDone,seed=seed,verbose=verbose)
#  mlefitHd = calcMLE(NOC,samples,popFreq,refData,condhd,knownRefhd, kit,modelDegrad,modelBWstutt,modelFWstutt,threshT,prC,lambda,fst,nDone=nDone,seed=seed,verbose=verbose)
#calcLRmcmc(mlefitHp,mlefitHd,niter = 20000,  quantile=0.1)

  #Show data:
  if(verbose) {
    if(isEPG) {
      plotEPG2(samples,kit,refData, AT=threshT) 
    } else {
      plotMPS2(samples,refData, AT=threshT) 
    }
  }
  #dat = prepareData(samples,refData,popFreq,threshT=threshT,minF=minF, normalize = normalize) #obtain data to use for analysis
  #rmp = calcRMPfst(dat,POIind,condInd,fst=fst ) #non-scaled values
  #LRupper = -sum(log10(rmp)) #maximum attainable LR (log10 scale)

  #Interpretations:
  #A: Model search:
  modelSelTable = NULL
  if(findOptimalModel) {
    
    #Running through all configurations (uses contLikMLE/validMLEmodel functions):
    modelDegrad = unique(c(FALSE,modelDegrad))
    modelBWstutt = unique(c(FALSE,modelBWstutt))
    modelFWstutt = unique(c(FALSE,modelFWstutt))
    minNOC = sum(condhp>0)
    #condOrder=condhd;knownRefPOI=POIind;knownRel=NULL;ibd=c(1,0,0);
    
    searchList = contLikSearch( NOC= minNOC:NOC, modelDegrad,modelBWstutt,modelFWstutt,samples=samples,popFreq=popFreq,refData=refData,condOrder=condhd,knownRefPOI=POIind, prC=prC,threshT=threshT,fst=fst,lambda=lambda,kit=kit,nDone=nDone,seed=seed,verbose=verbose,alpha=alpha)
    adjAIC = searchList$outtable[,2] #obtain crietion
    optimInd = which.max(adjAIC)[1] #get index of optimal model. Use simpler model if "Tie"
      #set to optimal model:
    mlehp = searchList$hpfitList[[optimInd]]
    mlehd = searchList$hdfitList[[optimInd]]
    modelSelTable = searchList$outtable #get model selection table
    
  } else { #otherwise use current model
    #if(verbose) print("Performing Maximum likelihood estimation (optimize)...")
    #DEG<-BWS<-FWS <- FALSE
    mlehp = calcMLE(NOC,samples,popFreq,refData,condhp,NULL,       kit,modelDegrad,modelBWstutt,modelFWstutt,threshT,prC,lambda,fst,nDone=nDone,seed=seed,verbose=verbose)
    mlehd = calcMLE(NOC,samples,popFreq,refData,condhd,knownRefhd, kit,modelDegrad,modelBWstutt,modelFWstutt,threshT,prC,lambda,fst,nDone=nDone,seed=seed,verbose=verbose)
  }
  
  #Obtain caclulated LR:
  LRmleObj = calcLRmle(mlehp,mlehd)
  LRupper = LRmleObj$log10LRupper #get upper LR
  LRmle = LRmleObj$log10LR
  LRmarker = LRmleObj$log10LRmarker
    
  #B: MLE VALIDATION:
  if(verbose) print("Performing MLE validations...")
  validhp = validMLEmodel(mlehp,plottitle = "Hp",alpha=alpha,createplot = verbose, verbose = verbose)
  validhd = validMLEmodel(mlehd,plottitle = "Hd",alpha=alpha,createplot = verbose, verbose = verbose)
  
  #C: Deconvolution:
  if(verbose) print("Performing Deconvolution...")
  DChp = deconvolve(mlehp)
  DChd = deconvolve(mlehd)
#  DChd$table2
  
  #D: Show expected PHs with data:
  if(verbose) {
    if(isEPG) {
      plotTopEPG2(mlehp,DChp) #show fitted PH under Hp
      plotTopEPG2(mlehd,DChd) #show fitted PH under Hd
    } else {
      plotTopMPS2(mlehp,DChp) #MLEobj=mlehp;DCobj=DChp;grpsymbol="_";locYmax=TRUE;options=NULL
      plotTopMPS2(mlehd,DChd) #show fitted PH under Hd
    }
  }
  
  #E: Perform MCMC based inference:
  mcmcobj = calcLRmcmc(mlehp,mlehd,niter = mcmcIter)
  LRmcmcBayes = mcmcobj$log10LRbayes
  LRcons = mcmcobj$log10LRcons
  
  #mcmchp = contLikMCMC(mlehp, delta=mcmcobj$delta,niter=mcmcIter)
  #validMCMC(mcmchp) #diagnostic of mcmc 
  
  #F: Perform Integral based inference:
  #mlefitHp=mlehp;mlefitHd=mlehd;
  LRintBayes = calcLRint(mlehp,mlehd,verbose = TRUE)
  
  #LRintBayes$calcHp
  #LRintBayes$calcHd
  LRint = LRintBayes$log10LR #obtain integration based LR
  
  #G: Calculate non-contributors (CAREFUL WITH NUMBER)
  tippets = calcTippet(POIind,mlehp,mlehd,nTippets) #MLE is default

  return( list(LRmle=LRmle,LRcons=LRcons,LRint=LRint,LRmarg=LRmcmcBayes,DChp=DChp$table2,DChd=DChd$table2,validhp=validhp,validhd=validhd,LRmarker = LRmarker,modelSelTable=modelSelTable, tippets=tippets) )
    
}
