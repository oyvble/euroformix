#' @title efm_createReport
#' @author Oyvind Bleka
#' @description Creating report from stored efm session
#' @details The function is a graphical layer for the functions in the package euroformix. See ?euroformix for more information.
#' @param mmTK An environment created with the efm function (can be loaded from project file)
#' @param type Indicate type of report to store (EVID or DC/DB)
#' @param save whether to save to text-file
#' @export


#load("C:/Users/oyvbl/Dropbox/Forensic/euroformix0/EFM4validation/Victor/proj.Rdata");type="EVID"
efm_createReport = function(mmTK, type="EVID",save=TRUE) { #function for storing MLE estimates of fitted models
  #generateEFMreport(type="EVID")
  #type = "EVID" , "START", "DC", "DB"
  sig=4 #number of significant levels
  colps="\t" #separator type  #h = list(action="EVID")
  myform = function(x) format(x,digits=sig) #helpfunction for formating numbers
  myround = function(x) round(x,sig)
  
  if(type=="START") type="EVID" #assume EVID if project just started
  set <- get(paste0("set",type),envir=mmTK) #get all setup-object 
  
  #If still not found when loading project
  if(is.null(set)) set <- get("setDC",envir=mmTK) #Assume as DC
  if(is.null(set)) set <- get("setDB",envir=mmTK) #Assume as DB  
  popKitInfo <- get("selPopKitName",envir=mmTK) #selected kit and population for popFreq
  
  printMLE <- function(mlefit,hyp) {
    mle <- cbind(mlefit$thetahat2,sqrt(diag(mlefit$thetaSigma2))) #standard deviation
    txt0 <- paste0("\n\n-------Estimates under ",hyp,"---------\n")
    txt1 <- paste0(c("Param.","MLE","Std.Err."),collapse=colps)
    for(i in 1:nrow(mle)) txt1 <- paste0(txt1,"\n",paste0( c(rownames(mle)[i],myform(mle[i,])),collapse=colps) )
    
    txt2 <- paste0("\n\nlogLik=",myround(mlefit$loglik))
    txt2 <- paste0(txt2,"\nadj.logLik=",myround(mlefit$loglik-length(mlefit$thetahat))) #adjust for number of params
    #txt2 <- paste0(txt2, "\nLik=", getSmallNumber(mlefit$loglik,sig))#format(exp(mlefit$loglik),digits=sig))
    txt <- paste0(txt0,txt1,txt2)
    return(txt)
  }
  
  sig=4 #number of significant levels
  colps="\t" #separator type  #h = list(action="EVID")
  myform = function(x) format(x,digits=sig) #helpfunction for formating numbers
  myround = function(x) round(x,sig)
  
  if(type=="START") type="EVID" #assume EVID if project just started
  set <- get(paste0("set",type),envir=mmTK) #get all setup-object 
  
  #If still not found when loading project
  if(is.null(set)) set <- get("setDC",envir=mmTK) #Assume as DC
  if(is.null(set)) set <- get("setDB",envir=mmTK) #Assume as DB  
  popKitInfo <- get("selPopKitName",envir=mmTK) #selected kit and population for popFreq
  
  printMLE <- function(mlefit,hyp) {
    mle <- cbind(mlefit$thetahat2,mlefit$thetaSE) #standard deviation
    txt0 <- paste0("\n\n-------Estimates under ",hyp,"---------\n")
    txt1 <- paste0(c("Param.","MLE","Std.Err."),collapse=colps)
    for(i in 1:nrow(mle)) txt1 <- paste0(txt1,"\n",paste0( c(rownames(mle)[i],myform(mle[i,])),collapse=colps) )
    
    txt2 <- paste0("\n\nlogLik=",myround(mlefit$loglik))
    txt2 <- paste0(txt2,"\nadj.logLik=",myround(mlefit$loglik-length(mlefit$thetahat))) #adjust for number of params
    txt2 <- paste0(txt2,"\nNumber of evals: ",mlefit$nEvals) #Number of evaluations
    txt2 <- paste0(txt2,"\nTime usage (sec): ",ceiling(mlefit$time)) #Number of evaluations
    #txt2 <- paste0(txt2, "\nLik=", getSmallNumber(mlefit$loglik,sig))#format(exp(mlefit$loglik),digits=sig))
    txt <- paste0(txt0,txt1,txt2)
    return(txt)
  }
  
  printSET <- function(model) { #print settings used
    txt <- paste0("\nEvidence(s)=",paste0(names(model$samples),collapse="/"))
    txt <- paste0(txt,"\nMarkers=",paste0(names(model$popFreq),collapse="/"))
    txt <- paste0(txt,"\n\n-------Model options-------")
    txt <- paste0(txt,"\nDetection threshold=",paste0(model$threshT,collapse="/"))
    txt <- paste0(txt,"\nFst-correction=",paste0(model$fst,collapse="/"))
    txt <- paste0(txt,"\nProbability of drop-in=",paste0(model$prC,collapse="/"))
    txt <- paste0(txt,"\nHyperparam lambda=",paste0(model$lambda,collapse="/"))
    txt <- paste0(txt,"\nDegradation:", ifelse(model$DEG,"YES","NO")) 
    txt <- paste0(txt,"\nBackward Stutter:",ifelse(model$BWS,"YES","NO")) 
    txt <- paste0(txt,"\nForward Stutter:",ifelse(model$FWS,"YES","NO"))
    txt <- paste0(txt,"\nBackward Stutter prop. prior=",paste0(deparse(eval(model$priorBWS)),collapse="") )
    txt <- paste0(txt,"\nForward Stutter prop. prior=",paste0(deparse(eval(model$priorFWS)),collapse="") )

    #Also adding additional settings:    
    txt <-  paste0(txt,"\nAdjusted fragmenth-length for Q-allele:",ifelse(model$adjQbp,"YES","NO")) #added in v4.0.0
    txt <-  paste0(txt,"\nRare allele frequency (minFreq):",set$minFreq)  #added in v3.0.0
    txt <-  paste0(txt,"\nNormalized after impute: ", ifelse(is.null(set$normalize) || set$normalize==1,"Yes","No") )  #added in v3.0.0 
    return(txt)
  }
  
  printMOD <- function(model,hyp) { #print refs
    txt <- paste0("\n\n-------Hypothesis ",hyp,"---------")
    txt <- paste0(txt,"\nNumber of contributors: ",model$nC) #Number of contributors
    txt <- paste0(txt,"\nKnown contributors: ",paste0(names(model$refData[[1]])[which(model$condOrder>0)],collapse="/")) #conditional references
    if(length(model$knownRef)) txt <- paste0(txt,"\nKnown non-contributors: ",paste0(names(model$refData[[1]])[model$knownRef],collapse="/")) #conditional references
    if(length(model$knownRel)) txt <- paste0(txt,"\nAssumed relationship: Last contributor is a ",names(model$ibd)[1]," to reference ",names(model$knownRel)) #Relationship
    return(txt)
  }
  
  printFREQ <- function(popFreq) { #print freqs
    locs = names(popFreq)
    txt <- paste0("\n\n-------Frequency data---------")
    for(loc in locs) {
      tmp <- paste0(names(popFreq[[loc]]),"=",popFreq[[loc]]) #get allele names with freqs
      txt <- paste0(txt,"\n",loc,": ",paste0(tmp,collapse="/"))
    }
    return(txt)
  }
  printDATA <- function() { #print data (evidence/refs)
    popFreq = set$popFreqQ
    samples = set$samples
    locs = names(popFreq)
    refData = set$refDataQ #assume structure [[loc]][[ref]] done by efm GUI
    refNames = NULL
    if(!is.null(refData)) refNames = unique(unlist(lapply(refData,names))) #get reference names

    txt <- paste0("\n\n-------Evaluating data---------\n")
    txt <- paste0(txt,"\t", paste0(names(samples),collapse=",")," | ")
    if(!is.null(refData)) txt <- paste0(txt, paste0(refNames,collapse=",")," | ")
    txt <- paste0(txt, "Freqs.\n")
    
    for(loc in locs) { #for each marker
      #loc=locs[1]
      txt <- paste0(txt,"\n",loc,"\n")
      alleles = names(popFreq[[loc]]) #get allele outcome
      if(!is.null(refData)) refLoc = refData[[loc]]
      for(allele in alleles) {
        line <- paste0(allele,"\t","| ")
        
        #INCLUDE PEAK HEIGHTS FOR EACH REPLICATE
        PHs = NULL
        for(sample in names(samples)) {
          #sample=names(samples)[1]
          av = samples[[sample]][[loc]]$adata
          ind = av==allele
          PH = "   " #empty peak height by defualt
          if(any(ind)) PH = samples[[sample]][[loc]]$hdata[ind] #get observed PH
          PHs = c(PHs,PH)
        }
        line = paste0(line,paste0(PHs,collapse=" ")) 
        
        #INCLUDE INDICATION OF WHETHER REFERENCE HAS ALLELE
        if(!is.null(refData)) {
          line <- paste0(line,"\t","| ")
          cross = NULL
          for(ref in refNames) {
            #ref=names(refLoc)[1]
            cross0 = " " #empty peak height by defualt
            if(ref%in%names(refLoc)) { #ensure that ref-name is included
              av = unlist(refLoc[[ref]])
              if( allele%in%av) cross0 = "x" #indicate with cross
            }
            cross = c(cross,cross0)
          }
          line = paste0(line,paste0(cross,collapse="   ")) 
        }
        
        #INCLUDE FREQUENCY LAST
        line = paste0(line,"\t |", signif(popFreq[[loc]][allele],5)) #small roundoff 
        txt <- paste0(txt,line,"\n") #add line
      } #end for each allele
    } #end each marker
    return(txt)
  }
  
  #START CURATING TEXT
  version = utils::packageVersion("euroformix")
  txt <- paste0("EuroForMix version ",version)#," (euroformix_",packageVersion("euroformix"),").")
  Rversion = strsplit(R.version.string," ")[[1]][3] #get only version number
  txt <- paste0(txt,"\nR-version: ",Rversion) #Add R-version used
  txt <- paste0(txt,"\nUser: ",Sys.getenv("USERNAME"))
  txt <- paste0(txt,"\nCreated: ",round(Sys.time()),"\n")
  
  txt <- paste0(txt,"\n-------Data-------")
  txt <- paste0(txt,"\nSelected STR Kit: ",popKitInfo[1])
  txt <- paste0(txt,"\nSelected Population: ",popKitInfo[2])
  
  txt <-  paste0(txt,printSET(set$mlefit_hd$model)) #Print Data and model options under Hd
  
  txt <- paste0(txt,"\n\n-------Optimalisation setting-------")
  txt <- paste0(txt,"\nRequired number of optimizations: ",set$mlefit_hd$nDone) #Added v3.0.0: Number of identical optimization
  txt <- paste0(txt,"\nAccuracy of optimisations (steptol): ",set$mlefit_hd$steptol) #Added v3.1.0: Steptol to use in nlm
  #txt <- paste0(txt,"\nSeed for optimisations: ", ifelse(is.null(set$mlefit_hd$seed),"NONE",set$mlefit_hd$seed)) #Added v3.0.0: 
  txt <- paste0(txt,"\nRestriction threshold (genotype outcome): ",set$mlefit_hd$resttol) #Added v4.1.0
  
  if(!is.null(set$mlefit_hp)) txt <- paste0(txt,printMOD(model=set$mlefit_hp$model,hyp="Hp")) #Print hypothesis Hp:
  if(!is.null(set$mlefit_hd)) txt <- paste0(txt,printMOD(model=set$mlefit_hd$model,hyp="Hd")) #Print hypothesis Hd:
  
  #store Hp,Hd-estimates 
  if(!is.null(set$mlefit_hp)) txt <- paste0(txt,printMLE(set$mlefit_hp$fit,"Hp"))
  if(!is.null(set$mlefit_hd)) txt <- paste0(txt,printMLE(set$mlefit_hd$fit,"Hd"))
  
  #store LR-estimates TO REPORT
  if(!is.null(set$mlefit_hp)) {
    res <- get(paste0("res", type),envir=mmTK) #extract correct result type resEVID/DB/DC
    
    #LR values:
    #txt1 <- paste0("MLE (sub-sub source): log10LR=",myform(log10(res$LRmle)))
    #txt2 <- paste0("MLE (sub source): log10LR=",myform(log10(res$adjLRmle))) 
    txt1 <- paste0("LR=",myform(res$LRmle))
    txt2 <- paste0("log10LR=",myform(log10(res$LRmle))) 
    txt3 <- paste0("Upper boundary: log10LR=",myform(log10(res$LRupper)))
    #txt4 <- paste0("log10LR (Laplace approximation)=",log10(res$LRlap)) 
    txt4 <- paste0(paste0(names(res$LRi),colps,colps,myform(res$LRi)),collapse="\n")
    txt <- paste0(txt,"\n\n-------MLE based LR (all markers)------\n",txt1,"\n",txt2,"\n",txt3,"\n")
    txt <- paste0(txt,"\n-------MLE based LR (per marker)------\n",txt4,"\n")
    
    #Obtain and insert NUMBER OF FAILED PP-plot points outside envelope: 
    alpha=get("optMLE",envir=mmTK)$alpha  #obtained from settings
    nFailedHp = res$nFailedHp
    nFailedHd = res$nFailedHd
    if(is.null(nFailedHp)) nFailedHp = NA #not computed
    if(is.null(nFailedHd)) nFailedHd = NA #not computed
    txtValid = paste0("Under H",c("p","d"),": ",c(nFailedHp,nFailedHd))
    
    if(!is.na(nFailedHp) || !is.na(nFailedHd)) {
      txt <- paste0(txt,"\n-------Model validation------\nNumber of fails (signif level=",alpha,"):\n",txtValid[1],"\n",txtValid[2])
    }
    
    #Print model searcher
    if(!is.null(set$SEARCH)) {
      txt1 <- paste0( colnames(set$SEARCH) ,collapse=colps) #obtain colnames
      for(i in 1:nrow(set$SEARCH)) txt1 <- paste0(txt1,"\n",paste0( set$SEARCH[i,] ,collapse=colps) )
      txt <-  paste0(txt,"\n\n-----Table of model comparisons-----\n",txt1) #add information 
    }
    
    #store consLR - estimate
    if(!is.null(res$MCMC)) {
      perc = res$MCMC$quantile*100 #obtain used percentile
      txt1 <- paste0("Conservative LR (",perc,"%): log10LR=",myform(res$MCMC$log10LRcons)) #insert conservative LR
      txt2 <- paste0("95% CI of conservative LR (",perc,"%): log10LR=[",paste0( myform(res$MCMC$log10LRconsCI),collapse=","),"]") #insert conservative LR
      txt3 <- paste0("Bayes Factor (MCMC): log10LR=",myform(res$MCMC$log10LRbayes)) #insert Bayesian LR
      txt4 <- paste0("Number of MCMC samples (setting): ",res$MCMC$nSamples)
      txt5 <- paste0("Variation of randomizer (setting): ",res$MCMC$delta) #selected delta
      txt6 <- paste0("Tuned variation of randomizer (estimated): ",signif( res$MCMC$deltaTuned,5) ) #Tuned delta
      txt7 <- paste0("Seed of randomizer (setting): ",res$MCMC$seed) 
      txt <- paste0(txt,"\n\n---RESULTS BASED ON MCMC SAMPLING---\n",txt1,"\n",txt2,"\n",txt3,"\n",txt4,"\n",txt5,"\n",txt6,"\n",txt7)
    }
    
    #store Bayes Factor - estimates
    if(!is.null(res$INT)) {
      txt1 <- paste0("Bayes Factor (integral): log10LR=",myform(log10(res$INT$LR))) #insert Bayesian LR
      err = myform(log10(res$INT$LRerror))
      txt2 <- paste0("Relative Error log10LR=[",err[1],",",err[2],"]") #insert Bayesian LR
      txt3 <- paste0("Relative Error (setting): ",res$INT$reltol)
      txt4 <- paste0("Deviation (setting): ",res$INT$dev) 
      txt5 <- paste0("MaxEvals (setting): ",res$INT$maxEval) 
      txt <- paste0(txt,"\n\n---RESULTS BASED ON NUMERICAL INTEGRATION---\n",txt1,"\n",txt2,"\n",txt3,"\n",txt4,"\n",txt5)
    }
    #store non-contribtor analysis
    if(!is.null(res$TIPPET)) {
      qqres = res$TIPPET$stat[[4]] #obtain quantile results
      
      txt1 = paste0("log10LR for model type: ",ifelse(res$TIPPET$type=="MLE","MLE approach","Bayesian approach") )
      txt2 = NULL
      for(i in seq_along(qqres)) txt2 = paste0(txt2,names(qqres)[i],"=",myround(qqres[i]),"\n")
      txt <-  paste0(txt,"\n\n-----RESULTS FROM NON-CONTRIBUTOR ANALYSIS-----\n",txt1,"\n",txt2) #add information 
      
    }
  }
  
  #Print data which are used for evaluation:
  txt <-  paste0(txt,printDATA()) #using info in set 
  
  #Print allele freqs last: ADDED in v2.0.1 (notice that only Hd is printed)
  #txt <-  paste0(txt,printFREQ(set$mlefit_hd$model$popFreq)) 
  #txt <-  paste0(txt,"\nRare allele frequency (minFreq):",set$minFreq)  #added in v3.0.0
  #txt <-  paste0(txt,"\nNormalized after impute: ", ifelse(is.null(set$normalize) || set$normalize==1,"Yes","No") )  #added in v3.0.0 
  
  txt <- as.matrix(txt)
  colnames(txt) <- paste0("This is a generated report from")
  if(save) tableSaver(txt , "txt") 
  return(txt)
} 
  