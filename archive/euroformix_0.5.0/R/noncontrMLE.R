#' @title noncontrMLE
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description noncontrMLE Calculate random log-likelhood values by exchanging a reference with a random man from population
#' @details The function simulates ntippet numbder of 
#' @param mlefit A fitted mle object from contLikMLE function
#' @param tipind Index in condOrder which are replaced with a random man
#' @param ntippet Number of samples to draw
#' @return TRUE/FALSE A boolean whether the likelihood vil be zero.
#' @export logL random log-likelihood values

noncontrMLE <- function(mlefit,tipind=NULL,ntippet=100) { 
  #tipref is index in refData to exchange with random man from population
     mod <- mlefit$model
     txt <- "Can't do tippet-analyis for given model"
     if(is.null(mod$refData) | is.null(mod$condOrder) | is.null(tipind))  stop(txt)
     if( is.na(mod$condOrder[tipind]) || mod$condOrder[tipind]==0) stop(txt)
     nU <- mod$nC - sum(mod$condOrder>0) #number of unknowns under Hp                    
     Glist <- getGlist(mod$popFreq) #get random man-Glist 
     refData <- mod$refData

     locs <- names(mod$refData) #loci to evaluate
     refind <- which(mod$condOrder>0) #conditional references under Hp
     refind <- refind[!refind%in%tipind] #remove tippet-ref  
    
     logL <- rep(-Inf,ntippet) #vector of log-likelihood values
     hpZero  <- rep(FALSE,ntippet) #boolean whether likelihood is zero under hp
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
          if(mod$prC==0 && any(hpZero[simind]==FALSE) ) hpZero[simind] <- hpZero[simind] | iszerolik(evid0,ref0,nU,mod$xi) #if no drop-in assumed
        }
       }
     }
     print(paste0("Optimizing ",sum(!hpZero)," likelihood values..."))
     for(m in 1:ntippet) { #for each random individual from the population
       if(!hpZero[m]) {
        for(loc in locs)  refData[[loc]][[tipind]] <-  Gsim[[loc]][m,]
        logL[m] <- contLikMLE(mod$nC,mod$samples,mod$popFreq,refData,mod$condOrder,mod$knownref,mod$xi,mod$prC,mlefit$nDone,mod$threshT,mod$fst,mod$lambda,delta=mlefit$delta,pXi=mod$pXi,kit=mod$kit,verbose=FALSE)$fit$loglik 
       }
       if(m%%(ntippet/10)==0) print(paste0(m/ntippet*100,"% finished..."))
     }
     return(logL)
} #end Tippet function

