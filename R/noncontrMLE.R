#' @title noncontrMLE
#' @author Oyvind Bleka
#' @description noncontrMLE calculate random log-likelhood values by exchanging a reference with a random man from population
#' @details The function simulates ntippet non-contributors
#' @param mlefit A fitted mle object from contLikMLE function (assuming Hp fitted object)
#' @param tipind Index in condOrder which are replaced with a random man
#' @param ntippet Number of samples to draw
#' @param verbose Boolean whether printing optimization progress. Default is FALSE.
#' @param seed The user can set seed if wanted
#' @param maxThreads Maximum number of threads to be executed by the parallelization

#' @return logL log-likelihood values of noncontributors
#' @export 

noncontrMLE <- function(mlefit,tipind=NULL,ntippet=100,verbose=FALSE,seed=NULL,maxThreads=32) { 
  if(!is.null(seed)) set.seed(seed) #set seed if provided
  
  #tipref is index in refData to exchange with random man from population
     mod <- mlefit$model
     txt <- "Couldn't do tippet-analyis for given model"
     if(is.null(mod$refData) | is.null(mod$condOrder) | is.null(tipind))  stop(txt)
     if( is.na(mod$condOrder[tipind]) || mod$condOrder[tipind]==0) stop(txt)
     nU <- mod$nC - sum(mod$condOrder>0) #number of unknowns under Hp                    
     Glist <- getGlist(mod$popFreq) #get random man-Glist 
     refData <- mod$refData

     locs <- names(mod$refData) #loci to evaluate
     refind <- which(mod$condOrder>0) #get conditional index of references under Hp
     refind <- refind[!refind%in%tipind] #remove tippet-ref  
    
     logL <- rep(-Inf,ntippet) #vector of log-likelihood values
     Gsim <- list()
     for(loc in locs) { #sample random individuals and check if they give Lik=0
       condR <- unlist(refData[[loc]][refind] ) #take out known refs under Hp 
       Gsim[[loc]] <-  Glist[[loc]]$G[ sample(1:length(Glist[[loc]]$Gprob),ntippet,prob=Glist[[loc]]$Gprob,replace=TRUE) ,] #Sample random genotypes from popFreqQ
       if(ntippet==1) unGsim <- t(Gsim[[loc]])
       if(ntippet>1) unGsim <- unique(Gsim[[loc]]) 
       for(j in 1:nrow(unGsim)) {
        ref0 <- c(unGsim[j,],condR) #conditional references
        simind <-  which(Gsim[[loc]][,1]==unGsim[j,1] & Gsim[[loc]][,2]==unGsim[j,2]) #get index of matching genotypes
        for(ss in names(mod$samples)) {
          evid0 <- mod$samples[[ss]][[loc]]$adata
        }
       }
     }
     print(paste0("Optimizing ",ntippet," likelihood values..."))
     for(m in 1:ntippet) { #for each random individual from the population
        for(loc in locs)  refData[[loc]][[tipind]] <-  Gsim[[loc]][m,] #fill in genotypes of reference
        logL[m] <- contLikMLE(nC=mod$nC,samples=mod$samples,popFreq=mod$popFreq,refData=refData,condOrder=mod$condOrder,knownRef=mod$knownref,
          xi=mod$xi,prC=mod$prC,threshT=mod$threshT,fst=mod$fst,lambda=mod$lambda,pXi=mod$pXi,kit=mod$kit,xiFW=mod$xiFW,pXiFW=mod$pXiFW,
          nDone=mlefit$nDone,delta=mlefit$delta,seed=mlefit$seed,verbose=FALSE)$fit$loglik 
         #arguments not used; maxIter=100,knownRel=NULL,ibd=c(1,0,0),seed=NULL,maxThreads=32
        
       if(m%%(ntippet/10)==0) print(paste0(m/ntippet*100,"% finished..."))
     }
     return(logL)
} #end Tippet function

