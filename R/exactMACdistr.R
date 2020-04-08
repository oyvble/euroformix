#' @title exactMACdistr
#' @author Oyvind Bleka
#' @description Calculation of false positive probability of Matching allele counter.
#' @details This function will calculate the false positive probability for a given range of allele matches.
#'
#' @param si A vector with each element as the sum of allele frequencies for observed alleles.
#' @return ret A vector with false positive probabilities for each given number of observed alleles.
#' @export


exactMACdistr <- function(si) {
 #si - A vector: With elements equal the sum of the observed frequencies in Evidence
 Wlist <- list()
 I <- length(si)
 for(i in 1:I) Wlist[[i]] <- cbind(0:2,c((1-si[i])^2,2*si[i]*(1-si[i]),si[i]^2))

 propagM <- numeric()
 for(i in 1:I) {
  mac <- Wlist[[i]][,1]
  pmac <- Wlist[[i]][,2]
  if(i>1) {
   newP <- outer(pmac,propagM[,2])   
   newM <- outer(mac,propagM[,1],FUN="+") #Changed: outer now includes function FUN="+" to avoid use of log and exp.
   agg <- aggregate(c(newP),by=list(c(newM)),sum)   #calculate marginals for updated vector again:
  } else {
   agg <- Wlist[[i]]
  }
  propagM <- cbind(agg[,1],agg[,2]) 
 }
 vec <- propagM[,2]
 names(vec) <- propagM[,1]
 return(vec)
}