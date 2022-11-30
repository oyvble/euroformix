#' @title prepareData
#' @author Oyvind Bleka 
#' @description Reorganasing data for calculation preparations
#' @details Optionally to use beforehand (still required for Qual. calcs)
#' Helpfunction to reorganising data which are input for further calculations (includes Qassignation). 
#' Also takes detection thresholds as argument to remove possible alleles below the thresholds.
#' 
#' @param mixData A list with evidence profiles [[sample]][[locus]]$adata/hdata
#' @param refData A list with reference profiles [[reference]][[locus]]$adata
#' @param popFreq A list with population frequencies [[locus]]
#' @param minF The freq value included for new alleles (new alleles as potential stutters will have 0). Default NULL is using min.observed in popFreq.
#' @param normalize Whether normalization should be applied or not. Default is FALSE.
#' @param threshT A detection threshold value or thresholds per makers (marker names must be defined). NULL ignores filtering.
#' @param fillHomGen Whether to fill in homozygote genotypes given with one allele
#' @param verbose Whether printing out information
#' @return Restructured data in list used as input for functions evaluating the liklihood function
#' @export

prepareData = function(mixData,refData=NULL,popFreq=NULL,minF=NULL,normalize=FALSE,threshT=NULL,fillHomGen=TRUE, verbose=TRUE) { #Helpfunction to get data to analyse
  if(is.null(popFreq)) stop("Populatation frequency object must be provided")
  locs <- names(popFreq) #get loci in popFreq (decides order)
  mixData2 <- lapply(mixData,function(x) return(x[locs])) #default is all data included
  if(!is.null(threshT)) { #if detection tresholds givven 
    for(loc in locs)  {
      AT = getMarkerVal(threshT,loc) #get potential marker specific threshold 
      
      for(mix in names(mixData)) { #for each evidence
        av = mixData2[[mix]][[loc]]$adata #get alleles
        hv = mixData2[[mix]][[loc]]$hdata #get heights
        if(length(av)>0) { #if contains alleles empty
          keep = hv>=AT
          av = av[keep]
          hv = hv[keep]
        } 
        mixData2[[mix]][[loc]] = list(adata=av,hdata=hv) #get updated vals
      }
    } #end for each loci
  } #if thresholds are provided

  refData2 = NULL
  if(!is.null(refData)) {
    refData2 <- list()
    for(loc in locs)  {
      refData2[[loc]] <- lapply(refData,function(x) {
        av = unlist(x[[loc]])#$adata #get alleles
        if(fillHomGen && length(av)==1) av = rep(av,2) #alleles given only ones m
        return(av) #return selected loci
      })
    }
  }
  
  Qret <- Qassignate(mixData2, popFreq, refData2,incS=FALSE,incR=FALSE, minF = minF, normalize = normalize,verbose=verbose)
  return(list(samples=Qret$samples,refData=Qret$refData,popFreq=Qret$popFreq))
}