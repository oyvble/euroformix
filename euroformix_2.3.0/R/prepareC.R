#' @title prepareC 
#' @author Oyvind Bleka
#' @description prepareC is used in the functions contLikMLE,contLikMarg and contLikMCMC to prepare input to C-call
#' @details The function builds the data input to the C-code
#' @param nC Number of contributors in model.
#' @param samples A List with samples which for each samples has locus-list elements with list elements adata and hdata. 'adata' is a qualitative (allele) data vector and 'hdata' is a quantitative (peak heights) data vector.
#' @param popFreq A list of allele frequencies for a given population.
#' @param refData Reference objects with list element [[s]][[i]]. The list element has reference-list with list-element 's' having a loci-list adata with list-element 'i storing qualitative data.
#' @param condOrder Specify conditioning references from refData (must be consistent order). For instance condOrder=(0,2,1,0) means that we restrict the model such that Ref2 and Ref3 are respectively conditioned as 2. contributor and 1. contributor in the model. 
#' @param knownRef Specify known non-contributing references from refData (indices). For instance knownRef=(1,2) means that reference 1 and 2 is known non-contributor in the hypothesis. This affectes coancestry correction.
#' @param kit shortname of kit: {"ESX17","ESI17","ESI17Fast","ESX17Fast","Y23","Identifiler","NGM","ESSPlex","ESSplexSE","NGMSElect","SGMPlus","ESX16", "Fusion","GlobalFiler","ForenSeq"}
#' @param knownRel Specify the index of the related contributing reference from refData (one index). For instance knownRel=2 means that unknown1 is related to reference 2 with ibd specified relationship.
#' @param ibd the identical by decent coefficients of the relationship (specifies the type of relationship)
#' @param fst The co-ancestry coefficient. Default is 0.
#' @param incS A boolean whether to include potential stutters to the allele outcome. Default is FALSE.
#' @return ret A list of data input to call the C-code with
#' @export 

prepareC = function(nC,samples,popFreq,refData=NULL,condOrder=NULL,knownRef=NULL,kit=NULL,knownRel=NULL,ibd=c(1,0,0),fst=0,incS=FALSE){
 Qallele="99" #Name of allele given if missing in evidence. Defualt is 99. This is important when considering the degradation model since 99 is closest to maximum allelein a locus. 
 LUSsymbol="_" #a character symbol used to separate repeatunit and LUS.
 #Note: Supports Invariant order of markers and caseletters!!
 #Supports replicates in the samples. 
 #should support to have zero contribution in some markers (all dropped out)-> all y=0
 #Evaluates only intercept of mixture-markers and freq-markers
 nS <- length(samples) #number of replicates
 #test 1) Require same markers in all replicates. 
 locs <- lapply(samples,function(x) toupper(names(x)))  #marker names for each replicates
 nLs <- sapply(locs,length)
 nL <- unique(nLs)
 if(length(nL)>1) stop("Number of markers in the replicates was not the same")
 locs <- unique(unlist(locs)) #get unique markernames
 if(nL!=length(locs)) stop("Different markers was specified in the replicates")
 names(popFreq) <- toupper(names(popFreq)) #toupper case!
 locs <- locs[locs%in%names(popFreq)] #take intercept with locus in evidence and population-freq
 #print( paste0("Evaluated loci: ", paste0(locs,collapse=",") ) )
 if(nC<1) stop("Number of contributors was specified as less than 1!")

 if(is.null(condOrder)) {
  condOrder <- rep(0,nC) #insert condorder if missing
  nK <- 0 #number of known contributors
 } else { #check that they are unique and in right order
  tmp <- condOrder[condOrder>0]
  nK <- length(tmp) #number of known contributors
  if( nK!=length(unique(tmp)) ) stop("Specify unique positions!")
  if( any( sort(tmp,decreasing=FALSE)!=(1:nK)) ) stop("Please condition references starting from 1. position")
 }
 
 #convertion of values in popFreq, mixData and Glist$G:
 #loci-order follows as in mixData: "locs". Rearrange names:
 popFreq <- popFreq[locs] #order popFreq to mixData-order
 
 #if(!is.null(knownRef) && !is.null(knownRel) && ibd[1]<1 && knownRel%in%knownRef) knownRef = setdiff(knownRef,knownRel) #activated in v2.3.0. THIS IS THE LRMIX MODEL: ENSURE THAT REF GIVEN AS RELATED IS NOT GIVEN AS KNOWN NON-CONTRIBUTOR.
 #else: warning("The typed related individual was also included as a known non-contributor. This means that he/she is typed twice! This will deviate from LRmix!")
 if( !is.null(knownRel) && ibd[1]<1 && fst>0 && !knownRel%in%c(knownRef,which(condOrder>0)) ) error("The related reference must be included as either a contributor or a known non-contributor! Please re-specify the hypothesis.")  #activated in v2.3.0: be sure that knownRel is already in knownRef (the non-contributor under Hd). Alternative  #knownRef = c(knownRef,knownRel) , but then also must considered as known non-contributor under Hp

 #Get probability of 1st unknown for all genotype outcomes (may contain related)
 Gset <- Gprob <- list() #store in a list for each markers
 for(loc in locs) {
  refK = unlist(refData[[loc]][knownRef]) #extract alleles of prev. typed individuals
  if(!is.null(condOrder)) refK = c(refK, unlist(refData[[loc]][which(condOrder>0)])) #extract alleles of prev. typed individuals
  refR = unlist(refData[[loc]][knownRel]) #extract alleles of related individual (must be length 2 if given!) 
  Glist <- calcGjoint(freq=popFreq[[loc]],nU=1,fst=fst,refK=refK,refR=refR,ibd=ibd)
  Gset[[loc]] <- Glist$G 
  Gprob[[loc]] <- Glist$Gprob
 }

 #Fix references: Assign condition to condM-matrix
 condM <- matrix(-1,nrow=nL,ncol=nC) #default is no references (=-1)
 #assign references to condM-matrix by values of Glist
 if(!is.null(refData) && !is.null(condOrder) && any(condOrder>0)) {
  names(refData) <- toupper(names(refData)) #toupper case!
  for(loc in locs) {
   subRef <- refData[[loc]] #take out relevant reference
   if(length(subRef)==0)  stop(paste('Missing locus (',loc,') in refData.',sep=''))
   for(k in 1:length(subRef)) { #for each reference
    if(length(subRef[[k]])==0) next #updated from v.0.6.4: Allowing unknown contributors for missing markers.
    if(length(subRef[[k]])!=2) stop("References need to have exactly two alleles.") #updated from v.0.6.4
    if(condOrder[k]>0) {
     Gind1 <- subRef[[k]][1]==Gset[[loc]][,1] & subRef[[k]][2]==Gset[[loc]][,2]
     Gind2 <- subRef[[k]][2]==Gset[[loc]][,1] & subRef[[k]][1]==Gset[[loc]][,2]
     condM[which(loc==locs),condOrder[k]] = which(Gind1 | Gind2) - 1 #subtract with one since we work from 0-indice
    }
   }
  }
 }

 #NEXT STEPS DO FOLLOWING:
 #Step 1) Add potential non-observed stutters to popFreq w/allele freq 0
 #Step 2) Create the vector BPind indicating what allele(index) the ith allele stutter to (given allele order in popFreq).
 #Step 3) Encode allelenames to 0:(nA-1) with nA=number of alleles
 isLUS0 <- any(sapply(samples,function(x) sapply(x,function(y) any(grepl(LUSsymbol,y$adata)) ))) #detect whether LUS variant is used. "_" needed in any of the alleles
 allLUS0 <- all(sapply(samples,function(x) sapply(x,function(y) all(grepl(LUSsymbol,y$adata)) ))) #check that all alleles are LUS format
 if(isLUS0 && !allLUS0) { #IN CASE OF ONLY SOME ALLELES WITH LUS VARIANT
   stop("Only some alleles contained the LUS syntax _. All or no alleles must contain this syntax. Program stops!")
 }

 #HELPFUNCTION USED FOR LUS:
 getLUSstutter = function(x) { #get the stuttered variants of LUS
   tmp = as.numeric(x) 
   stutt2 <- paste0(tmp[1]-1,LUSsymbol,tmp[2]-1) #stutter-allele
   if(length(tmp)>2) stutt2 <- paste0(c(stutt2,tmp[-(1:2)]),collapse=LUSsymbol) #add other LUS variants
   return(stutt2)
 }

 popFreq2 <- popFreq #new version with decoded alleles
 #UPDATE: ADD POTENTIAL STUTTERS TO popFreq with frequency=0. This is first time where potential stutters should be included (can also be included in Qassignate to get compatible results for previous EFM versions).
 #decode allele-names to index names + Stutter-preparation
 #Potential stutters first be added to both popFreq (old) and popFreq2 (new), since genotypes must be recognized 
 allASind <- as.numeric() #vector with allele indices (starts from index 1) used to denote what allele an allele backward-stutters to.
 for(loc in locs) { #for each marker
#loc=locs[23]
   anames <- names(popFreq[[loc]]) #old names
   ASind <- rep(0,length(anames)) #default is no stutter indicing 
   if(!all(anames%in%Qallele) && incS) { #If no system dropout (i.e. only allele 99 is present) or stutters are included
    #STEP 1: Add potential non-observed stutters to popFreq w/allele freq 0
    anames2 <- anames[anames!=Qallele]
    if( isLUS0 ) {    #Important: possible LUS variant must be handled at each step.
     splitA <- strsplit(anames2,LUSsymbol) #split the alleles
     stutt <-  sapply(splitA, getLUSstutter) #get stutter variants
    } else { #if it is a number
     stutt <- as.character(as.numeric(anames2)-1) #get stutter variant, convert back to character
    }
    stuttAadd <- stutt[!stutt%in%anames] #stutter alleles to add
    anames <- c(anames,stuttAadd)
    popFreq[[loc]] <- c(popFreq[[loc]],rep(0,length(stuttAadd))) #add zero freqs
    names(popFreq[[loc]]) <- anames

    #Identify what allele each allele will backwardstutter to (given by index)
    if( isLUS0 ) {
     indUse <- anames!=Qallele
     anames2 <- anames[indUse] #remove the Qallele
     stutt = rep(NA,length(anames))
     stutt[indUse] <- sapply( strsplit(anames2,LUSsymbol), getLUSstutter ) #get stuttering alleles
    } else {
     anames <- as.numeric(anames) #convert from string to numbers
     stutt <- anames-1 #get corresponding stutters
    } #end in case of LUS
    ASind <- match(as.character(stutt),as.character(anames),nomatch=0) #get index of what allele it stutters to (anames -> stutt)
   } #end if include Stutter
   allASind <- c(allASind,ASind) #add to list
  
   anames2 <- 0:(length(popFreq[[loc]])-1) #new names
   popFreq2[[loc]] <- popFreq[[loc]]
   names(popFreq2[[loc]]) <- anames2 #update names in popFreq
   #go through each observed alleles in samples and update names
   for(s in 1:nS) { #for each sample
    obsA <- samples[[s]][[loc]]$adata  #observed alllees
    if(length(obsA)>0) { #if atleast 1 observed
     for(j in 1:length(obsA)) { #for each observed allele we encode allele to indices
      if(!any(anames==obsA[j])) stop(paste0("For locus ",loc,": Please add allele ",obsA[j]," to the population frequncy table"))
      samples[[s]][[loc]]$adata[j] <- anames2[anames==obsA[j]] #update allele-name
     }
    } #dont do anything if non-observed!
   } #end for each samples
 } #end for each marker
 #fix genotypes:

 #Encode allel-names in Gset: NB: BE CAREFUL USE TEMP-VARIABLES!
 Gset2 <- Gset #keep an old version!
 for(loc in locs) { #for each marker
  oldnames<-names(popFreq[[loc]])
  newnames<-names(popFreq2[[loc]])
  for(old in oldnames) Gset[[loc]][old==Gset2[[loc]]] <- newnames[which(old==oldnames)] #change values
 }

 #Counting alleles in known/typed references: MUST BE DONE AFTER INSERTING POTENTIAL STUTTERS IN POPFREQ
 mkvec <- numeric() #number of times each alleles are sampled
 nkval <- numeric() #total number of sampled alleles in each marker
 for(loc in locs) {
  tmp <- rep(0, length(popFreq[[loc]]))
  if(!is.null(condOrder) & !is.null(refData)) {
   for(k in 1:length(condOrder)) {
    if(condOrder[k]!=0) { 
     ind <- which( names(popFreq[[loc]])%in%refData[[loc]][[k]] )
     tmp[ind] = tmp[ind] + (length(ind)==1) + 1 #add twice sampled if homozygote #Bug discovered and updated in v1.11.3
    }
   }
  }
  if(!is.null(knownRef)) { #known non-contributors
   for(k in knownRef) {
    ind <- which( names(popFreq[[loc]])%in%refData[[loc]][[k]] )
    tmp[ind] = tmp[ind] + (length(ind)==1) + 1   #add twice sampled if homozygote  #Bug discovered and updated in v1.11.3
   }
  }
  nkval <- c(nkval,sum(tmp)) #number of sampled (for each loci)
  mkvec <- c(mkvec,tmp) 
 }

 #take into account for replicates here!
 #count and vectorize alleles in mixtures here:
 nA <- obsA <- obsY <- numeric()
 for(loc in locs) { #for each marker
  for(s in 1:nS) { #for each replicates
   nA <- c(nA, length(samples[[s]][[loc]]$adata)) #count number of alleles 
   obsA <- c(obsA, samples[[s]][[loc]]$adata) #add a-observations
   obsY <- c(obsY, samples[[s]][[loc]]$hdata) #add h-observations
  }
 }
 CnA <- c(0,cumsum(nA)) #cumulative number of alleles for each marker

 #LAST: get kit-info and find size of alleles
 slist <- list() #size list used for degeneration
 mlist <- list() #mean of fragment size for each markers
 if(!is.null(kit)) { 
   kitinfo <- getKit(kit)
   if(length(kitinfo)>1) { #Updated: check that it is not only NA returned
    for(loc in locs) { #for each markers
     slist[[loc]] <- numeric()
     subkit <- subset(kitinfo,toupper(kitinfo$Marker)==loc) #control on dyer
     mlist[[loc]] <- mean(subkit$Size)
     Avec <- names(popFreq[[loc]])
     for(an in Avec) {  #for each allele
      an2 <- an
      isLUS <- grepl(LUSsymbol,an2) #detect whether LUS variant is used
      if( isLUS ) an2 <- as.numeric(unlist(strsplit(an ,LUSsymbol)))[1] #convert to numeric
      ind <- which(subkit$Allele==an2)
      suppressWarnings( {
        if(length(ind)==0) ind <- which.min(abs(as.numeric(subkit$Allele) - as.numeric(an2))) #nearest neighbour (in allele name)-> maximum of the alleles
      })
      if(length(ind)==0) ind <- length(subkit$Allele) 
#       if(length(ind)==0) ind <- which.min( abs(subkit$Size - mlist[[loc]]) ) #Input with mean of the sizes
      size <- subkit$Size[ind]
      slist[[loc]] <- c(slist[[loc]],size) #include size
     }
    } #end for each locs
    if( any(sapply(slist,length)==0) ) stop("Wrong kit specified! It didn't contain all markers given in data.")
   } else {
    stop("Wrong kit name specified! It was not recognized by getKit().")
   }
 }

 #dependent on popFreq and ref-conditions
 nG <- sapply(Gprob,length) #number of genotype combinations
 CnG <- c(0,cumsum(nG))
 CnG2 <- c(0,cumsum(nG*2)) #note: 2 columns for each genotype!!
 bp <- (unlist(slist)-125)/100 #vectorize allele base pair size over all loci, diff with minimum average bs for each marker and scale with 100 always
 pG <- unlist(Gprob) #vectorize genotype probabilities over all loci
 Gvec <- as.integer(rbind(unlist(Gset))) #vectorize a big matrix (loci are put chronologic)
 condRef <- c(condM) #vectorized over all loci
 nAall <- sapply(popFreq,length) #Number of population-alleles on each loci
 CnAall <- c(0,cumsum(nAall)) #cumulative number of alleles
 pA <- unlist(popFreq) #need each allele probability for drop-in probabilities

 retlist <- list(nC=as.integer(nC),nK=as.integer(nK),nL=as.integer(nL),nA=as.integer(nA), obsY=as.numeric(obsY),obsA=as.integer(obsA),CnA=as.integer(CnA),allASind=as.integer(allASind),nAall=as.integer(nAall),CnAall=as.integer(CnAall),Gvec=as.integer(Gvec),nG=as.integer(nG),CnG=as.integer(CnG),CnG2=as.integer(CnG2),pG=as.numeric(pG),pA=as.numeric(pA), condRef=as.integer(condRef),mkvec=as.integer(mkvec),nkval=as.integer(nkval),nS=as.integer(nS),bp=as.numeric(bp))
 return(retlist)
} #end function

