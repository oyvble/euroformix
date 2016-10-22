#' @title Qassignate
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description Q-designation. It also takes care of include new alleles into popFreq by assigning it lowest observed frequnce.
#' @details Assignes non-shown alleles as one single allele.
#' @param samples A List with samples which for each samples has locus-list elements with list elements adata and hdata. 'adata' is a qualitative (allele) data vector and 'hdata' is a quantitative (peak heights) data vector.
#' @param popFreq A list of allele frequencies for a given population.
#' @param refData Reference objects with list element [[s]]$adata[[i]]. The list element has reference-list with list-element 's' having a loci-list adata with list-element 'i storing qualitative data.
#' @param doQ A boolean whether to Q-designate or not
#' @param incR A boolean whether to include reference-alleles in the Q-assignation 
#' @param incS A boolean whether to include potential stutters in the Q-assignation 
#' @return ret A list(popFreq,refData) with Q-designated alleles. 
#' @export

Qassignate <- function(samples,popFreq,refData=NULL,doQ=TRUE,incS=TRUE,incR=TRUE) {
 locs <- names(popFreq)
 minF <- min(unlist(popFreq)) #lowest observed frequency
 for(loc in locs) { #make Q-assignation for each loci
  evid0 <- evid <- unique(unlist( lapply(samples,function(x) x[[loc]]$adata) )) #vectorize alleles for all replicates
  if(!is.null(refData) && incR) evid  <- unique(c(evid,unlist(refData[[loc]])))
  if(incS) evid  <- unique(c(evid,as.character(as.numeric(evid)-1)))
  #tmp <- c(as.character(as.numeric(evid0)+1),as.character(as.numeric(evid0)+2)) #not necessary?
  #if(incS) evid <-  unique(c(evid,tmp[tmp%in%names(popFreq[[loc]])]))
  #if new alleles in evidence or references not in popFreq
  newa <- evid[!evid%in%names(popFreq[[loc]])]   
  if(length(newa)>0) {
   tmp <- names(popFreq[[loc]])
   popFreq[[loc]] <- c(popFreq[[loc]],rep(minF,length(newa)))
   names(popFreq[[loc]]) <-  c(tmp,newa)
   print(paste0("Locus ",loc,": Allele(s) ",paste0(newa,collapse=",")," was inserted with frequency ",minF))
  }
  popFreq[[loc]] <- popFreq[[loc]]/sum(popFreq[[loc]]) #normalize

  if(doQ) {#if Q-assignate
   tmp <- popFreq[[loc]][names(popFreq[[loc]])%in%evid] #find observed alleles
   if(length(tmp)<length(popFreq[[loc]])) { 
    tmp <- c(tmp,1-sum(tmp))
    names(tmp)[length(tmp)] <- "99"
   }
   popFreq[[loc]] <- tmp
   if(!incR && !is.null(refData)) { #insert 99 as default allele of missing refs
    newP <- names(popFreq[[loc]]) 
    if(!all(unlist(refData[[loc]])%in%newP)) { #there was some missing alleles
     for(k in 1:length(refData[[loc]])) refData[[loc]][[k]][!refData[[loc]][[k]]%in%newP] <- "99" #insert missing     
    }
   }
  }
 }
 return(list(popFreq=popFreq,refData=refData))
}
