#' @title calcQual
#' @description Implementation of the qualitative likelihood calculations
#' @details Assumes that prD is a identical for all contributors. Otherwise running forensim::likEvid (or possibly other implementations)
#' From likEvid documentation: If there are two replicates, showing alleles 12,13, and 14 respectively, then evids should be given as c(12,13,0,14), where the 0 is used as a separator. An empty replicate is simply 0. For example, replicates (12,13) and and one empty replicate must be given as: c(12,14,0,0).
#' @param evids vector of alleles present at a given locus for any number of replicates. Same as 'Repliste' in likEvid function. 
#' @param knownContr vector of genotypes for the known contributors. Genotype 12/17 should be given as a vector c(12,17) and genotypes 12/17,14/16, should be given as a unique vector: c(12,17,14,16). If empty, set to 0.
#' @param knownNonContr vector of genotypes for the known non-contributors
#' @param nUnknowns Number of unknown individuals under H. Set to 0 if there are no unknown contributors.
#' @param theta theta correction, value must be taken in [0,1)
#' @param prD probability of heterozygotes dropout. It is possible to assign different values per contributor. 
#' @param prC probability of drop-in applied per locus
#' @param freq vector of the corresponding allele frequencies of the analysed locus in the target population
#' @return likelihood value
#' @export 

calcQual = function(evids, knownContr, knownNonContr,nUnknowns,theta, prD, prC,freq) {
  
  #Check if conventional forensim must be run
  calcWithForensim = FALSE #whether forensim should be used directly
  if( length(prD) > 1) {
    prDunique = unique(prD)
    if( length(prDunique)>1 ) { 
      calcWithForensim = TRUE #need to use forensim if different values
    } else {
      prD = prDunique #can use own implementation if same value
    }
    #ALSO MUST CALC IF PrDhom is not pDhet^2
  } #if else( length(prDHom)>1 || prDHom != prDHet^2) 
    
  if(calcWithForensim) {
    nContr = length(knownContr)/2 + nUnknowns
    pDvec = rep(prD,nContr)
    #lik = forensim::likEvid(evids,knownContr,knownNonContr,nUnknowns,theta, pDvec,pDvec^2,prC,freq)
    warning("Not implemented!")
    lik = NA
    return(lik)    
  }
  
  #CONTINUE WITH OWN IMPLEMENTATION
  alleleNames = names(freq) #obtain allele names
  nAlleles = length(freq) #number of alles
  
  #CREATE A BOOLEAN OF WHETHER EVIDENCE HAS ALLELE IN FREQUENCY DATA
  hasEvidMat = NULL
  hasEvidVec0 = rep(0L,nAlleles)
  replicateIdx = 1
  alleleCounter=0 #number of visited alleles
  for(i in seq_along(evids)) {
    allele = evids[i]
    alleleCounter = alleleCounter + 1 #update
    if(alleleCounter==1) hasEvidVec = hasEvidVec0 #reset vector
    if(allele>0) { #if observed allele
      hasEvidVec[allele==alleleNames] = 1L #indicate as observd
    } else if(alleleCounter>1) { #otherwise if more than 1 allele is visited
      replicateIdx =  replicateIdx + 1 #skip to next replicate
      alleleCounter = 0
      hasEvidMat = cbind(hasEvidMat,hasEvidVec)
    }
  }
  hasEvidMat = cbind(hasEvidMat,hasEvidVec) #always append
  hasEvidVec = as.integer(hasEvidMat)
  nReplicates = ncol(hasEvidMat) #number of replictes
  
  #Prepared typed info:
  nTyped <- contrAlleles <- rep(0,nAlleles) #number of typed alleles (Each type)
  for(i in seq_len(nAlleles)) {
    if(!is.null(knownContr)) {
      numAlleles = sum(knownContr==alleleNames[i])
      if(numAlleles>0) {
        nTyped[i] = nTyped[i] + numAlleles
        contrAlleles[i] = contrAlleles[i] + numAlleles
      }
    }
    if(!is.null(knownNonContr)) nTyped[i] = nTyped[i] + sum(knownNonContr==alleleNames[i])
  }
  
  #Only running C++ function if at least 1 unknown:
  if(nUnknowns>0) {
    lik = .C("calcQual1marker", as.numeric(0), as.numeric(prD),as.numeric(prC),as.numeric(theta), as.integer(nUnknowns), as.integer(nAlleles), as.integer(nReplicates), as.integer(contrAlleles), as.integer(nTyped),  as.numeric(freq), as.integer(hasEvidVec))[[1]]
  } else {
    
    #NO GENOTYPE TRAVERSION NEEDED WHEN nUnknowns=0
    hasContr = contrAlleles > 0 #bool of which has contributin
    hasNoContr = !hasContr
    dropAlleleContr = prD^contrAlleles#[hasContr] #calculate log  prD^n_a (used two places)
    
    likEvid = 1
    for(r in seq_len(nReplicates)) {
      hasEvid = hasEvidMat[,r] #obtain observed indication for replicate
      A_contrSet = which(hasEvid & hasContr) #Contribution set
      B_dropoutSet = which(!hasEvid & hasContr) #Dropout set
      C_dropinSet = which(hasEvid & hasNoContr) #Dropin set
      if(length(A_contrSet)>0) likEvid = likEvid*prod(1-dropAlleleContr[A_contrSet])
      if(length(B_dropoutSet)>0) likEvid = likEvid*prod(dropAlleleContr[B_dropoutSet])
      if(length(C_dropinSet)>0) {
        likEvid = likEvid*prod(prC*freq[C_dropinSet])
      } else {
        likEvid = likEvid*(1-prC) #no dropin
      }
    }
    lik = likEvid
  }
  return(lik)
} #end function

