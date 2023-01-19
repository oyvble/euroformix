#' @title calcGjoint
#' @author Oyvind Bleka
#' @description getGlist Returns a list of joint genotypes with corresponding joing probabilities for unknowns contributors for a given marker.
#' @details The function returns the list of all possible genotypes, with corresponding probabilities. The allele-names in popFreq needs to be numbers.
#' If a related unknown is provided that will always be the first one
#' @param freq A vector of allele frequencies for a given population.  
#' @param nU Number of unknowns (where the LAST is a related)
#' @param fst Assumed theta/fst-correction
#' @param refK contains a vector of alleles for the known typed reference profiles (a,b,c,d...)
#' @param refR contains a vector of alleles for a related reference profile (a,b)
#' @param ibd the identical by decent coefficients of the relationship denotation
#' @return Glist A list with genotypes and genotype probabilities 
#' @export 
#' @examples
#' freq = rgamma(8,1,1)
#' freq = stats::setNames(freq/sum(freq),1:length(freq))
#' foo1 = calcGjoint(freq,nU=2,fst=0.1,refK=c("2","3","1","1"))
#' foo2 = calcGjoint(freq,nU=2,fst=0.1,refK=c("2","3","1","1"),
#'  refR=c("2","3"),ibd=c(1/4,1/2,1/4))

calcGjoint = function(freq,nU=1,fst=0,refK=NULL,refR=NULL,ibd=c(1,0,0)) {
 if(length(fst)!=1) stop("Wrong input length for fst")
 if(length(nU)!=1) stop("Wrong input length for number of unknowns (nU)")
 if(nU==0) stop("You must specify at least one unknown to use this function!")
 if(!is.numeric(freq))  stop("freq argument must be numeric!") 
 
 sumsToOne = all.equal(1,sum(freq)) #checking if freqs sums to one (this must always be the case for EFM)
 if(is.character(sumsToOne)) warning("freq argument must sum to one!")
 
 #Function calculates genotypes for all possib
 nn = length(freq) #number of alleles
 nG <- nn*(1+nn)/2 #number of allele outcome
 #nG^nU #NUMBER OF ITERATIONS
 #print(nG^nU)

 #PRESTEP: GET GENOTYPE OUTCOME: SIMILAR TO getGlist function
 av <- names(freq)   
 suppressWarnings({   
   if(!any(is.na(as.numeric(av))) )  av <-  as.numeric(av) #convert to numbers if 
 })

 #Modified in v2.3.0: Gmatrix is the vectorized upper triangular (1,1),(1,2),...,(1,n),(2,2),(2,3),...,(2,n),....,(n,n)
 G = numeric()
 for(i in 1:nn) {
  G = rbind(G, cbind( av[rep(i,nn - i + 1)], av[i:nn] ))
 }
 if(nrow(G)!=nG) warning("The size of the genotype outcome is not as expected!") #Must be true!
 ishomG <-  G[,1]==G[,2] #find G variants which are homozygous
 
 #COUNT NUMBER OF ALLELES IN refK:
 mkvec <- rep(0,length(av)) #count number of typed alleles
 if(!is.null(refK) && length(refK)>0 ) {
  tab <- table(unlist(refK))
  mkvec[match(names(tab),av)] <- tab  #add counted alleles
  names(mkvec) <- av
 }
 
 #####################################
 #GENERAL FORMULA WITH KINSHIP MODULE#
 #####################################
 
 P = function(Pi,ni,n) (ni*fst + (1-fst)*Pi)/(1+(n-1)*fst) #helpfunction for allele prob
 genoProb = function(mkvec2,refR=NULL,ibd2=NULL) { #helpfunction to get genotype prob for spec
   Gprob1 <- Gprob0 <- rep(NA,nG) #used to store genotype prob|typed alleles (all outcome)
   ishomR <- refR[1]==refR[2] #is the ref homozygous
   
   #prepare Gprobs for K0:
   for(i in 1:nG) { #for each genotype
     refG <- G[i,] #get proposed genotype
     indf <- match(refG,names(freq)) #get (allele) index of frequency
     tmpval = P(freq[indf[1]],ni=mkvec2[indf[1]],n=sum(mkvec2)) #prob. of 1st allele
     if(ishomG[i]) { 
       Gprob0[i] <- tmpval*P(freq[indf[2]],ni=mkvec2[indf[2]]+1,n=sum(mkvec2)+1) #get hom freq
     } else {
       Gprob0[i] <- 2*tmpval*P(freq[indf[2]],ni=mkvec2[indf[2]],n=sum(mkvec2)+1) #get het freq   
     }
   } #end for each type
   #sum(Gprob0) #MUST BE 1
   if(is.null(refR)) return(Gprob0) #return if no relatedness
   
   #Get frequencies
   A1ind <- match(G[,1],av) #index of allele1 prob, all(freq[A1ind]==pA1 )
   A2ind <- match(G[,2],av) #index of allele2 prob, all(freq[A2ind]==pA2 )
   nA <- as.integer(av%in%refR) #get typed alles of reference refR
   if(ishomR) nA[nA==1] <- 2 #should be 2 if homozygous
   
   #Get allele sharing types (between G and refR):
   ind1 <- G[,1]%in%refR
   ind2 <- G[,2]%in%refR
   mac <- as.numeric(ind1)+as.numeric(ind2)
   if(!ishomR) mac[ishomG & mac==2] <- 3 #set special case when Ghom has both allele in Rref when R is het variant
   
   #FOR NON-SHARING VARIANTS (K0 only)
   ind <- mac==0
   Gprob1[ind] <- Gprob0[ind]*ibd2[1]  #must multiply by ibd
   
   #FOR SHARING ONE ALLELE VARIANT (K0+K1)
   ind <- mac==1
   indUse <- rep(NA,sum(ind)) #This is index of alleles NOT IN REF (for situations one sharing allele)
   indUse[ind2[ind]] <-  A1ind[ind & ind2] 
   indUse[ind1[ind]] <-  A2ind[ind & ind1] 
   
   if(ishomR) {
     Gprob1[ind] <- Gprob0[ind]*ibd2[1]  + P(freq[indUse],ni=mkvec2[indUse],n=sum(mkvec2))*ibd2[2]
   } else {
     Gprob1[ind] <- Gprob0[ind]*ibd2[1]  + P(freq[indUse],ni=mkvec2[indUse],n=sum(mkvec2))/2*ibd2[2] 
   }
   ind <- mac==3 #ONE SHARING BUT WITH G AS homozygous variants 
   Gprob1[ind] <- Gprob0[ind]*ibd2[1]  + P(freq[A1ind[ind]],ni=mkvec2[A1ind[ind]],n=sum(mkvec2))/2*ibd2[2]  #USING ONLY pA1 always OK?
   
   #FOR SHARING BOTH ALLELE VARIANTS (K0+K1+K2 only)
   ind <- mac==2 #this will be only one of the variants
   if(ishomR) {
     Gprob1[ind] <- Gprob0[ind]*ibd2[1] + P(freq[A1ind[ind]],ni=mkvec2[A1ind[ind]],n=sum(mkvec2))*ibd2[2] + ibd2[3] #formula ok for ref hom/het
   } else {
     Gprob1[ind] <- Gprob0[ind]*ibd2[1] + (P(freq[A1ind[ind]],ni=mkvec2[A1ind[ind]],n=sum(mkvec2)) + P(freq[A2ind[ind]],ni=mkvec2[A2ind[ind]],n=sum(mkvec2)))/2*ibd2[2] + ibd2[3] #formula ok for ref hom/het
   }
   #sum(Gprob1) #MUST BE 1
   
   return(Gprob1)
 }  
 #CREATE RECURSION FUNCTION TO CALCULATE pG for a given counted alleles and prev. given pG
 Gprob <- array(data = NA, dim = rep(nG,nU)) #create structure. Each index correspond to genotype in G
 
 #function returning genoprob for given start mkvec and pGeno
 calcGprob = function(mkvecTMP,track=NULL,Uk=1,pGenoTMP=1) { #Uk gives depth in recursion. Stops after finishing loop on Uk=1 (i.e. returning from all)
   #track=NULL for Uk=1
   
   #EVALUATE RELATED AS LAST CONTRIBUTOR
   refR2 <- ibd2 <- NULL #default is no relatedness
   if(Uk==nU) { #if last contributor
     refR2 <- refR
     ibd2  <- ibd
   }
   pGenotmp <- pGenoTMP*genoProb(mkvecTMP,refR2,ibd2) #update probabilities conditioned on prev typed alleles
   
   for(gind in seq_len(nG) )  { #for each genotype
     mkvectmp <- mkvecTMP #obtain a new copy here
     
     for(l in 1:2) { #UPDATE ALLELE COUNTER:
       indf <- which(G[gind,l]==av) #get index of frequency	 
       mkvectmp[indf] <-  mkvectmp[indf] + 1  #update counter of sampled allele
     }
     
     if(Uk<nU) { #if still more unknowns to calculate for
       calcGprob(mkvecTMP=mkvectmp,track=c(track,gind),Uk=Uk+1,pGenoTMP=pGenotmp[gind] )  #recurse deeper
     } else { #final unknown has been achievied: Need to store genotype and return back
       #print(paste0(paste0(track,collapse="-"),"-",gind,":",pGenotmp[gind])) #print tracks
       Gprob[ rbind(c(track,gind)) ] <<- pGenotmp[gind] #store calculations
     } 
   } #end for each genotype
 }
 #RUN RECURSION
 calcGprob(mkvecTMP=mkvec)
 #sum(Gprob) #PROB OF G OUTCOME MUST BE 1  

 return(list(G=G,Gprob=Gprob)) #return list
}


