#' @title freqs_listToTable
#' @author Oyvind Bleka
#' @description Converting frequency data from list to table format
#' @param popFreq List with frequencies per marker
#' @return Table with alleles ordered. Ready to be exported to frequency file. 
#' @export
#' 
freqs_listToTable = function(popFreq) {
  nL <- length(popFreq)
  unAchr <- unique(unlist(lapply( popFreq,names) )) #get character alleles
  ord <- order(as.numeric(unAchr))  #all alleles must be able to convert to numeric
  unAchr <- unAchr[ord]  #make increasing order
  
  outtab = matrix("",ncol=nL,nrow=length(unAchr)) #table to save to file
  for(i in 1:nL) { #go through all markers
    freqs <- popFreq[[i]] #get frequencies
    outtab[ match(names(freqs),unAchr) ,i ] = freqs #insert frequencies
  } 
  outtab = cbind(unAchr,outtab)
  colnames(outtab) = c("Allele",names(popFreq)) #insert marker names
  return(outtab)  
} #end of functions

