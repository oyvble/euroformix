#' @title getGlist 
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description getGlist Returns a list of genotypes with corresponding probabilities for a given allele frequency-list.
#' @details The function returns the list of all possible observed genotypes, with corresponding Hardy-Weinberg assumed probabilities. The allele-names in popFreq needs to be numbers.
#' @param popFreq A list of allele frequencies for a given population. Each named element in the list must be a allele-named vector with allele frequencies. 
#' @return Glist A list with genotypes and genotype probabilities for each locus.
#' @export 

 getGlist <- function(popFreq) {
  locs <- names(popFreq)
  Glist <- list()
  for (loc in locs) {
   G = t(as.matrix(expand.grid(rep(list(as.numeric(names(popFreq[[loc]])),as.numeric(names(popFreq[[loc]])))))))
   keep = G[2, ] >= G[1, ]
   G <- G[, keep]
   G <- matrix(as.character(G), nrow = 2)
   tmpP = t(as.matrix(expand.grid(rep(list(as.numeric(popFreq[[loc]]),as.numeric(popFreq[[loc]]))))))
   tmpP <- tmpP[, keep]
   if(is.null(dim(tmpP))) tmpP <- cbind(tmpP)
   Gprob = exp(colSums(log(tmpP)))
   ishet = G[1, ] != G[2, ]
   Gprob[ishet] = 2 * Gprob[ishet]
   Glist[[loc]] <- list(G = t(G), Gprob = Gprob)
  }
  return(Glist)
 }
