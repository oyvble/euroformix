#' @title simDOdistr
#' @author Oyvind Bleka
#' @description MCMC Allele dropout distribution sampler based on total number of alleles in an evidence.
#' @details simDOdistr samples from the drop-out distribution based on total number of alleles in evidence under a specified prepositions. It returns if no samples was accepted in first iteration
#' @param totA Total number of allele-observations in evidence.
#' @param nC Number of contributors to assume in the preposition.
#' @param popFreq Population frequencies listed for each loci popFreq[[locname]]
#' @param refData List with alleles to condition on. ref[[locname]][[referencename]]=(A1,A2)
#' @param M The number of samples for each iteration.
#' @param minS The number of minimum accepted samples.
#' @param prC Assumed drop-in probability. Can be a vector (must contain the marker names)
#' @param pDknown A vector of known drop-out probabilities for each contributors. Default is NA which means it is unknown.
#' @return Vector with accepted samples from the dropout distribution
#' @export
#' @examples
#' \dontrun{
#' popFreq <- list() #create population frequencies
#' for(i in 1:3) {
#'  freqs <- rgamma(rpois(1,10),1,1)
#'  popFreq[[paste0("loc",i)]] <- stats::setNames(freqs/sum(freqs),1:length(freqs))
#' }
#' simDOdistr(6,2,popFreq)
#' }

simDOdistr= function(totA,nC,popFreq,refData=NULL,M=1e3,minS=2000,prC=0,pDknown=rep(NA,nC)) {
 #totA - total number of alleles observed
 #nC - number of contributors (in each hypotheses)
 #popFreq - population frequencies for each loci
 #prC - drop-in probability
  
 primtall = .getPrim() #get prime numbers
 
 #Prepare data:
 refvec <- list()
 locs <- names(popFreq) #loci to consider
 nL <- length(locs)
 
 #Check dropin prob and prepare dropin list
 if(length(prC)==1) {
   pCv = rep(prC,nL) #common dropin prob for all markers
 } else {
   pCv = setVecRightOrder(prC,  locs) 
 }
 prC_list = list() #store drop-in list
 for(i in 1:nL) {
   pos = ceiling(log(1e-16)/log(pCv[i])) #maximum number of position
   prC_vec =  1-pCv[i]/(1-pCv[i])
   if(pos>0) prC_vec <- c(prC_vec,pCv[i]^c(1:pos)) #dropin probabilities
   prC_list[[i]] = prC_vec
 }
 
 uH <- rep(NA,nL) #number of unknowns
 for(loc in locs) {
  i <- which(locs==loc)
  refvec[[loc]] <- numeric()
  uH[i] <- nC
  if(!is.null(refData)) {
   refs <- refData[[loc]]
   for(j in 1:length(refs)) {
    refvec[[loc]] <- c(refvec[[loc]], refs[[j]])
    uH[i] <- uH[i] - as.integer(length(refs[[j]])>0) #number of unrestricted unknowns: Updated in v1.9  
   }
  }
 }
 if(any(uH<0)) { print("There was more references than contributors"); return(NULL) }

 #Uknown numbers under hyps
 done<-FALSE
 cPrD_dist <- NULL

 while(!done) {
  aCount = rep(0,M)
  prDvec = runif(M) #sample prDs
  tmpcount = rep(0,M)
  for(loc in locs) {
   i <- which(locs==loc)
   prC_vec = prC_list[[i]] #get list of dropin
   freq <- popFreq[[loc]]
   freqN = names(freq) #updated: not convert to numeric

   HRef = t(matrix(rep(refvec[[loc]],M),ncol=M))
   HA2 <- cbind(HRef,matrix(sample(freqN,2*M*uH[i],freq,replace=TRUE),nrow=M))
   HA <- matrix(1,ncol=ncol(HA2),nrow=nrow(HA2))
   for(j in 1:length(freqN))  HA[HA2==freqN[j]] = primtall[j]    #rename with primenumbers

   #generate dropout:
   Z_h = matrix(runif(M*ncol(HA)),ncol=ncol(HA)) #get Z uniform samples
   Z_h2 = matrix(FALSE,nrow=M,ncol=ncol(HA)) #default: all are false
 
   if(all(is.na(pDknown))) { #can treat every column with same drop-out probabiltiy
    Z_h2[Z_h < prDvec] <- TRUE #TRUE if dropped out
   } else { #need to go through each double-columns seperately
    for(j in 1:(ncol(Z_h)/2)) {
     colind <- (2*j-1) : (2*j)
     if(is.na(pDknown[j])) {
      Z_h2[,colind] <- Z_h[,colind] < prDvec
     } else if(pDknown[j]>0){ #drop-out given
      Z_h2[,colind] <- Z_h[,colind] < pDknown[j]
     }
    }
   }
 
   #Let droppouts get primenumber=1. These are always equal something
   for(j in 1:ncol(HA))  HA[ Z_h2[,j] ,j] = 1

   #idea: Knows that product of different primes is unique
   #takes product of previous unique visited primers. Dropouts are also checked for!
   HA_tmp = HA[,1]
   cc_H = as.numeric(HA_tmp!=1) # counter: Don't count if dropout(equal 1)
   for(j in 2:ncol(HA)) { 
    equal = HA_tmp%%HA[,j]==0 #| HA[,j]==1 #check if next primer is in product (not a new unique number)
    cc_H[!equal] = cc_H[!equal] + 1 #add to counter if not equal
    HA_tmp[!equal] = HA_tmp[!equal]*HA[!equal,j] #get product of primes
   }
  
   #generate dropin:
   #idea:
   # 1) generate number of dropouts
   # 2) For each size-number of contamination we generate prime numbers and
   #check the number of the prime-products to see if the cont-alleles are unique

   #samples number of contaminations in each case:
   ncontam = sample(0:(length(prC_vec)-1),M,prob=prC_vec,replace=TRUE) #sample number of dropin
   contmax = max(ncontam)
   if(contmax>0) { #if any contanimation
    for(j in 1:contmax) { #for each size of contamination:
      sel = ncontam>=j
      if(sum(sel)>0) { #total number of contaminations
       prim = primtall[sample(1:length(freq),sum(sel),prob=freq,replace=TRUE)] #sample kind
       equal = HA_tmp[sel]%%prim==0
       cc_H[sel][!equal] = cc_H[sel][!equal] + 1 #add to counter if not equal
       HA_tmp[sel][!equal] = HA_tmp[sel][!equal]*prim[!equal] #get product of primes
     }
    } #end for each cont. size
   } #end if cont

   #count up uniques for loci 
   tmpcount = tmpcount  + cc_H 
  } #end for:iind #each loci

  cPrD_dist <- c(cPrD_dist,prDvec[tmpcount==totA])
  if(length(cPrD_dist)==0) return(cPrD_dist) #return if no samples was found in first iteration
  if(length(cPrD_dist)>=minS) done = TRUE
 } #end while
 return( cPrD_dist ) 
}
