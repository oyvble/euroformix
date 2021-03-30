#' @title Qassignate
#' @author Oyvind Bleka
#' @description Q-designation. It also takes care of include new alleles into popFreq by assigning it lowest observed frequnce.
#' @details Assignes non-shown alleles as one single allele.
#' @param samples A List with samples which for each samples has locus-list elements with list elements adata and hdata. 'adata' is a qualitative (allele) data vector and 'hdata' is a quantitative (peak heights) data vector.
#' @param popFreq A list of allele frequencies for a given population.
#' @param refData Reference objects with list element [[s]]$adata[[i]]. The list element has reference-list with list-element 's' having a loci-list adata with list-element 'i storing qualitative data.
#' @param doQ A boolean whether to Q-designate or not
#' @param incR A boolean whether to include reference-alleles in the Q-assignation 
#' @param incS A boolean whether to include potential BW stutters in the Q-assignation 
#' @param minF The freq value included for new alleles (new alleles as potential stutters will have 0). Default NULL is using min.observed in popFreq.
#' @param normalize A boolean of whether normalization should be applied or not. Default is FALSE.
#' @param verbose A boolean of whether printing out information
#' @return ret A list(popFreq,refData,samples) with Q-designated alleles. 
#' @export

Qassignate <- function(samples,popFreq,refData=NULL,doQ=TRUE,incR=FALSE,incS=FALSE,minF=NULL,normalize=FALSE,verbose=TRUE) {
 Qallele="99" #Name of allele given if missing in evidence. Default is 99. This is important when considering the degradation model since 99 is closest to maximum allelein a locus. 
 LUSsymbol="_" #a character symbol used to separate repeatunit and LUS.
 locs <- names(popFreq)
 if(is.null(minF)) minF <- min(unlist(popFreq)) #lowest observed frequency if given as NULL
 for(loc in locs) { #make Q-assignation for each loci
  evid <- unique(unlist( lapply(samples,function(x) x[[loc]]$adata) )) #vectorize alleles for all replicates

  #if new alleles in evidence or references: insert minFreq and normalize
  if(!is.null(refData) && incR) evid  <- unique(c(evid,unlist(refData[[loc]]))) #update evid to also include ref-alleles

   #Add n-1 stutters to the evid vector. BLOCK keeped to be compatible with previous EFM versions (before v3)
   if(length(evid)>0 && incS) { #added: Must have observed allele to consider stutter
    isLUS <- grepl(LUSsymbol,evid) #detect whether LUS variant is used
    evidS <- numeric()
    if( all(isLUS) ) {
      tmp <- t( matrix(as.numeric(unlist(strsplit(evid,LUSsymbol))),ncol=length(evid)) )#convert to numeric
      stutt <- paste0(tmp[,1]-1,LUSsymbol,tmp[,2]-1)
      if(ncol(tmp)>2) { #in case of additional LUS-variants: The last variants are added
       stutt2 <- apply(tmp[,-(1:2),drop=F],1,function(x) paste0(x,collapse=LUSsymbol) )
       stutt <- paste0(stutt,LUSsymbol,stutt2) 
      }
      evid <- unique(c(evid,stutt))
    } else if(!any(isLUS)) { #case of no LUS
     evid  <- unique(c(evid,as.character(as.numeric(evid)-1)))
    }  else {
    stop(paste0("Only some of the alleles in locus ",loc," contained the LUSsymbol=",LUSsymbol))
   }
  }

  newa <- evid[!evid%in%names(popFreq[[loc]])]   #get alleles not in popFreq-table
  
  if(length(newa)>0) {
   tmp <- names(popFreq[[loc]])
   popFreq[[loc]] <- c(popFreq[[loc]],rep(as.numeric(minF),length(newa)))
   names(popFreq[[loc]]) <-  c(tmp,newa)
   if(verbose) print(paste0("Locus ",loc,": Allele(s) ",paste0(newa,collapse=",")," was inserted with frequency ",minF))

   if(as.logical(normalize)) { #Update in v2.0: Normalization is now an option
    popFreq[[loc]] <- popFreq[[loc]]/sum(popFreq[[loc]]) #normalize
   }
  }

  if(doQ) {#if Q-assignate, i.e. setting non-observed alleles of references as Qallele ("99")
   tmp <- popFreq[[loc]][names(popFreq[[loc]])%in%evid] #find observed alleles
   
   #Assign a Q-allele only if frequency table is greater than those obtained in evidence
   if( length(tmp) < length(popFreq[[loc]])) {  
    tmp <- c(tmp,1-sum(tmp))
    names(tmp)[length(tmp)] <- Qallele
   }
   popFreq[[loc]] <- tmp
   if(!incR && !is.null(refData)) { #insert 99 as default allele of missing refs
    newP <- names(popFreq[[loc]])  #new population table (alleles)
    if(!all(unlist(refData[[loc]])%in%newP)) { #there was some missing ref-alleles in population table
      
      #V3.3.0: Introduce Q-allele to population table if not already defined (using rare freq assigning)
      if(!Qallele%in%newP) { #only if Q-allele not defined
        popFreq[[loc]] = setNames( c(popFreq[[loc]],as.numeric(minF)), c(names(popFreq[[loc]]),Qallele)  )
        if(verbose) print(paste0("Locus ",loc,": Q-allele was inserted with frequency ",minF))
        if(as.logical(normalize)) popFreq[[loc]] <- popFreq[[loc]]/sum(popFreq[[loc]]) #normalize
      }
      
      for(k in 1:length(refData[[loc]])) refData[[loc]][[k]][!refData[[loc]][[k]]%in%newP] <- Qallele #insert missing     
    }
   }
  }
 }
 return(list(popFreq=popFreq,refData=refData,samples=samples))
}
