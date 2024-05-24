#' @title prepareC
#' @author Oyvind Bleka
#' @description Used for preparing C++ calls
#' @param nC Number of contributors in model.
#' @param samples A List with samples which for each samples has locus-list elements with list elements adata and hdata. 'adata' is a qualitative (allele) data vector and 'hdata' is a quantitative (peak heights) data vector.
#' @param popFreq A list of allele frequencies for a given population.
#' @param refData Reference objects has locus-list element [[i]] with a list element 'r' which contains a 2 long vector with alleles for each references.
#' @param condOrder Specify conditioning references from refData (must be consistent order). For instance condOrder=(0,2,1,0) means that we restrict the model such that Ref2 and Ref3 are respectively conditioned as 2. contributor and 1. contributor in the model. 
#' @param knownRef Specify known non-contributing references from refData (index). For instance knownRef=(1,2) means that reference 1 and 2 is known non-contributor in the hypothesis. This affectes coancestry correction.
#' @param kit shortname of kit: Obtained from getKit()
#' @param BWS Boolean of whether back-stutter model should be used
#' @param FWS Boolean of whether for-stutter model should be used
#' @param AT The analytical threshold given. Used when considering probability of allele drop-outs.
#' @param pC A numeric for allele drop-in probability. Default is 0.
#' @param lambda Parameter in modeled peak height shifted exponential model. Default is 0.
#' @param fst is the coancestry coeffecient. Default is 0.
#' @param knownRel gives the index of the reference which the 1st unknown is related to.
#' @param ibd the identical by decent coefficients of the relationship (specifies the type of relationship)
#' @param minF The freq value included for new alleles (new alleles as potential stutters will have 0). Default NULL is using min.observed in popFreq.
#' @param normalize Whether normalization should be applied or not. Default is FALSE.
#' @param adjFragQallele Indicate whether fragmenth length of Q-allele is based on averaged weighted with frequencies
#' @return ret A list of data input to call the C++ code with
#' @export 

prepareC = function(nC,samples,popFreq, refData, condOrder, knownRef, kit,BWS,FWS,
                 AT,pC,lambda,fst,knownRel,ibd, minF,normalize, adjFragQallele) {
  Qallele = "99" #substitution name for Q-alleels (non observed alleles)
  LUSsymbol="_" #a character symbol used to separate repeatunit and LUS.
  MPSsymbol = ":" #Added in version 3.1.0. Used for extracting CE for MPS strings. Example is "10:[ATCG]10".
  if(is.null(minF)) minF = min(unlist(popFreq))
  
  #HELPFUNCTION USED FOR LUS:
  getLUSstutter = function(x,isBW=TRUE) { #get the BW stuttered variants of LUS
    tmp = as.numeric(x) 
    stuttSign = -1 #BW is default
    if(!isBW) stuttSign = +1 #FW stutter
    stutt2 <- paste0(tmp[1] + stuttSign,LUSsymbol,tmp[2] + stuttSign) #BW stutter-allele
    if(length(tmp)>2) stutt2 <- paste0(c(stutt2,tmp[-(1:2)]),collapse=LUSsymbol) #add other LUS variants
    return(stutt2)
  }
  
  #Obtain replicate and marker names etc
  nReps = length(samples) #Number of samples overall
  repNames = names(samples) #obtain replicate names
  for(rep in repNames) names(samples[[rep]]) =  toupper(names(samples[[rep]])) #ensure UPPER CASE marker names
  locs_evids = lapply(samples,function(x) names(x)) #obtain locus names per marker (UPPER CASE)
  locs = unique(unlist(locs_evids)) #obtain unique markers
  locs = intersect(locs,names(popFreq)) #Consider loci only given in popFreq
  nLocs = length(locs)
  
  #Handle Reference data: (they may be in EFM format)
  isEFMformat = FALSE
  MAXnRefs = 0 #obtain maximum found references (across all markers)
  refNamesCond <- refNamesAll <- NULL #references to condition on
  if(!is.null(refData)) {
    if( any(toupper(names(refData))%in%locs) ) isEFMformat = TRUE #There were marker names in refData
    
    if(!isEFMformat) { #converting if not in EFM format:
      for(k in seq_along(refData)) names(refData[[k]]) = toupper(names(refData[[k]])) #ensure upper letter
      refData2 = list()
      for(loc in locs) refData2[[loc]] = lapply(refData,function(x) x[[loc]])
      refData = refData2
    } 
    names(refData) = toupper(names(refData)) #ensure upper locus name letters
    
    #Check references and standardize the name (refNames)
    nRefsLoc = sapply(refData,function(x) length(x)) #obtain maximum
    refNamesAll = names(refData[[which.max(nRefsLoc)]]) #obtain names
    refNamesUnique =  unique(unlist(lapply(refData,function(x) names(x)))) #obtain unique references
    if( !all(refNamesUnique%in%refNamesAll) ) stop("Some markers had missing reference profiles!")
    MAXnRefs = max(nRefsLoc)
  }  
  if(!is.null(knownRef)) if(any(knownRef>MAXnRefs)) stop("Indicated known non-contributor reference not found!")
  if(length(knownRel)>1) stop("Multiple related not supported!")
  
  if(!is.null(condOrder)) {
    if(any(condOrder>MAXnRefs) ) stop("Indicated contributors not found!")
    refNamesCond = refNamesAll[condOrder] #obtain reference names to condition on
  }

  #Indicate number of stutter models
  nStutterModels = 0
  if(BWS) {
    nStutterModels = 1
    if(FWS) nStutterModels = 2
  } else if(FWS) {
    stop("Model cannot assume only Forward Stutter Model (without assuming Backward Stutter Model)!")
  }
  nStutterModels = rep(nStutterModels,nLocs) #must be marker specific
  
  #Prepare hyp params
  if(length(AT)==1) {
    ATv = rep(AT,nLocs) #LOD
  } else {
    ATv = setVecRightOrder(AT,locs)
  }
  if(length(pC)==1) {
    pCv = rep(pC,nLocs) #dropout prob
  } else {
    pCv = setVecRightOrder(pC,locs)
  }
  if(length(fst)==1) {
    fstv = rep(fst,nLocs) #theta correction
  } else {
    fstv = setVecRightOrder(fst,locs)
  }
  if(length(lambda)==1) {
    lambdav = rep(lambda,nLocs) #lambda used for dropin
  } else {
    lambdav = setVecRightOrder(lambda,locs)
  }
  
  #OBTAIN Hypotheses information:
  NOK = sum(condOrder>0) #obtain number of references to condition on 
  knownGind <- matrix(-1,ncol=nLocs,nrow=NOK) #default is no references (=-1)
  relGind <- rep(-1,nLocs) #matrix(-1,ncol=nLocs,nrow=nC) #default is no references (=-1)
  ibd0 = c(1,0,0) #no relation
  if(!is.null(ibd)) ibd0 = ibd #assign
    
  #PREPARE DATA VECTORS:
  nAlleles <- nPotStutters <- nRepMarkers <- rep(0,nLocs) #num alleles and extra potential stutters  
  QalleleIndex <- rep(-1,nLocs) #default is no dropout
  YvecLong <- FvecLong <- DvecLong <- numeric()
  BWtoLong <- FWtoLong <- BWfromLong <- FWfromLong <- numeric() #backward/forward stutter info (to-allele)
  useDEG = FALSE #note: overrides input (useful for show topEPGplot)
  if(!is.null(kit)) { #always turn DEG on if kit name is provided
    kitinfo = euroformix::getKit(kit) #get kit information
    if(!is.null(kitinfo) || !is.na(kitinfo[1])) useDEG = TRUE #use degradation model
  }
  basepairLong <- numeric() #base pair information
  nTypedLong <- maTypedLong <- numeric() #= rep(0,nLocs)  #number of typed (total)
  alleleNames <- alleleNamesALL <- NULL #store name of all alleles
  
  nJointGenos <- nGenos <- nUnknowns <- nKnowns <- rep(0L,nLocs) #number ofgenotypes (also joint) to traverse
  genoList = list() #store alleles in genotypes
  for(m in seq_len(nLocs)) { #m=2
    loc = locs[m] #extract locus
    freqAll = popFreq[[loc]] #get freqs
    
    #Look on observations
    #First obtain replicates to consider:
    repNamesMarker = NULL
    for(rep in repNames) {
      if(loc%in%locs_evids[[rep]]) repNamesMarker = c(repNamesMarker, rep)
    } 
    nRepMarkers[m] = length(repNamesMarker) #number of replicates
    
    #Then check the observed Peaks for each replicate (check AT)
    rep_alleles = lapply(samples[repNamesMarker], function(x) x[[loc]]$adata[x[[loc]]$hdata>=ATv[m]])
    alleles = unique(unlist(rep_alleles)) #get alleles
    nAlleles[m] = length(alleles) #store number of unique alleles
    
    #Check if alleles are in LUS format (i.e. contains the _ separation symbol):   
    isLUS <- FALSE
    allelesNotQ = setdiff(alleles,Qallele) #remove alleles with Q-alleles (may be with old input)
    isSTRING = FALSE #Need to skip stutter model if alleles are not numbers (string)
    if(length(allelesNotQ)>0) {
      tmp = grepl(LUSsymbol,allelesNotQ)
      if(all(tmp)) { #if all alleles was with LUS syntax
        isLUS = TRUE  #LUS IS CONSIDERED FOR THE MARKER!
      } else if(any(tmp)) { #if only some alleles are in LUS
        stop(paste0("In marker ",loc,": Only some alleles contained the LUS syntax. All or no alleles must contain this syntax for a marker. Program stops!"))

      } else {
        #CHECK if alleles can be converted to numbers (if not we need to skip stutter model)
        suppressWarnings({
          check = as.numeric(allelesNotQ) #convert to numbers
          if(any(is.na(check))) isSTRING = TRUE #stop("Data format not supported for Stutter model!")
        })        
      }
    }
    if(isSTRING) nStutterModels[m] = 0 #set number of stutter models to zero

    #Creating a matrix for the peaks per replicate given per row
    yv = matrix(0,ncol= nAlleles[m],nrow=nRepMarkers[m],dimnames = list(repNamesMarker,alleles)) #create PH-matrix
    for(rep in repNamesMarker) {
      locDat = samples[[rep]][[loc]]
      indUse = locDat$hdata>=ATv[m] #check those to use
      if(any(indUse)) yv[rep, match(locDat$adata[indUse],alleles)] = locDat$hdata[indUse] #insert PHs
    }
    
    #Obtain allele frequences
    freqs = as.numeric()
    if( nAlleles[m] > 0) { #IF ANY ALLELES TO use (NOT EMPTY LOCI)
      newA = alleles[!alleles%in%names(freqAll)] #new alleles
      if(length(newA)>0) {
        newA = setNames(rep(minF,length(newA)),newA)
        freqAll = c(freqAll,newA)
      }
      if(normalize) freqAll = freqAll/sum(freqAll)
      freqs = freqAll[match(alleles,names(freqAll))]
    }
    freqsNotObserved = freqAll[!freqAll%in%names(freqs)] #obtain obtained non-observed alleles
    
    #CHECKING Q-allele
    Qfreq = 1-sum(freqs)
    names(Qfreq) = Qallele
    addQ = Qfreq>0 #Q-allele freq must be at last 0
    
    dropinWeight = numeric() #init empty
    if( nAlleles[m]>0 ) { #IF ANY ALLELES TO KEEP (NOT EMPTY LOCI)
      FvecLong = c(FvecLong , freqs) #add freqs for obsalleles
      freqs2 = t(replicate(nRepMarkers[m],freqs)) #span out for each replicate
      if(length(freqs)==1) freqs2 = t(freqs2) #need to rotate if only one allele in frequency
      dropinWeight = log(freqs2) + log(pCv[m]) + dexp(yv-ATv[m],rate=lambdav[m],log=TRUE)
      dropinWeight[is.infinite(dropinWeight)] = -1e+100 ##insert large negative number
    } 
    if(addQ) {
      QalleleIndex[m] = nAlleles[m] #last index is the Q-allele (adjusted for C++)
      FvecLong = c(FvecLong , Qfreq) #add freq for dropout
      
      #Add zero to last column (for Q-alleles)
      insZeros = rep(0,nRepMarkers[m]) #insert number of zeros
      yv = cbind(yv , insZeros) #add  zero PH as last column (Q-alleles)
      dropinWeight = cbind(dropinWeight,insZeros) #add  zero dropinweigt as last column (Q-alleles)
      nAlleles[m] = nAlleles[m] + 1 #number of alleles to traverse in genotypes (observed + dropout)
    } 
    DvecLong = c(DvecLong , as.numeric(dropinWeight)) #Add drop-in weights
    YvecLong = c(YvecLong , as.numeric(yv)) #add PHs for obs alleles
    
    #Obtain fragment sizes
    if(useDEG) { #Extracting basepairs
      kittab = kitinfo[toupper(kitinfo$Marker)==loc,,drop=FALSE] #obtain kit info for specific loci
      CEallelesAll = names(freqAll) #allelesNotQ #is already CE allele by default
      if(isLUS) {
        CEallelesAll <- as.numeric( sapply(strsplit(CEallelesAll ,LUSsymbol),function(x) x[1] ) ) #extracting first allele  
      } else if( all(grepl(MPSsymbol,CEallelesAll)) ) {
        CEallelesAll <- as.numeric( sapply(strsplit(CEallelesAll ,MPSsymbol),function(x) x[1] ) ) #extracting first allele  
      } 
      
      #Obtain bp for all alleles (also defined in freq data used for Q-allele)
      fragLengthAll = .getFragLength(CEallelesAll,kittab,isSTRING ) #obtain corresponding fragment length of all alleles (also freq)
      fragLengthObservedNotQ = fragLengthAll[match(allelesNotQ,names(freqAll))]
      basepairLong = c(basepairLong, fragLengthObservedNotQ) #Extend vector
      if(addQ) {
        QalleleFraglength = max(kittab$Size) #use longest fragment by default
        if(adjFragQallele) { #use weighted average
          allelesNotObserved = setdiff(names(freqAll),allelesNotQ)
          if(length(allelesNotObserved)>0) { #require at least one observation
            weight = freqAll[allelesNotObserved]
            weight = weight/sum(weight) #normalize
            QalleleFraglength = sum(fragLengthAll[allelesNotObserved]*weight) #obtain weighted expectation
          }
        }
        basepairLong = c(basepairLong,QalleleFraglength)
      }
    } 
    
    #Include Q-allele if not already added (and should be added)
    if(addQ) alleles = c(alleles,Qallele) #this must contain also Q-allele
    
    #OBtain matrix of alleles for each genotype outcome
    Gset = numeric()
    for(i in 1:length(alleles)) Gset = rbind(Gset, cbind( alleles[rep(i,length(alleles) - i + 1)], alleles[i:length(alleles)] ))
    genoList[[loc]] = Gset #add genotypelist

    #Prepare BWstutter/FWstutter relations:
    #isNum = !any(is.na(av)) #check if allele is num
    #isNum = isNum && (BWS || FWS) #also
    nPotStutters[m] = 0  #number of potential stutters (Alleles not part of av). Default is none
    BWalleles <- FWalleles <- NULL
    if( (BWS || FWS) && !isSTRING) {
      if(isLUS) { #CHECK FIRST IF LUS
        tmplist = strsplit(allelesNotQ ,LUSsymbol) #split up elements with LUS separator
        if(BWS) BWalleles <- sapply(tmplist, getLUSstutter,isBW=TRUE) #get corresponding BW stutters
        if(FWS) FWalleles <- sapply(tmplist, getLUSstutter,isBW=FALSE) #get corresponding FW stutters
      } else { #Otherwise it is CE based
        if(BWS) BWalleles = as.character(as.numeric(allelesNotQ)-1)
        if(FWS) FWalleles = as.character(as.numeric(allelesNotQ)+1)
      }
      PSalleles = setdiff(c(BWalleles,FWalleles),alleles)
      alleles2 = c(alleles,PSalleles) #include potientials stutters in additional to all
      nPotStutters[m] = length(PSalleles) #get number of potential stutters (all dropout)
      
      BWto = match( BWalleles,alleles2,nomatch=0 )-1 #get indices where each allele gives BW contribution to
      FWto = match( FWalleles,alleles2,nomatch=0 )-1 #get indices where each allele gives FW contribution to
      if(!FWS) FWto = rep(-1,length(BWalleles)) #not used
      BWfrom = match(as.character(alleles2),BWalleles,nomatch=0)-1 #get index of what allele it receive stutter from  (Allele 1=index0)
      FWfrom = match(as.character(alleles2),FWalleles,nomatch=0)-1 #get index of what allele it receive stutter from  (Allele 1=index0)
      
    } else { #ELSE IF STUTTER NOT USED
      BWto <- FWto <- rep(-1,length(allelesNotQ)) 
      BWfrom <- FWfrom <- rep(-1,length(alleles))  
      alleles2 = alleles #copy allele names
    }
    if(addQ) {  #last index is dummy variable (Dropout doesn't stutter). THerefor need QalleleIndex to not be used
      BWto = c(BWto, -1)
      FWto = c(FWto, -1) 
    }
    alleleNames = c(alleleNames, alleles) #store alleles (including potential stutters
    alleleNamesALL = c(alleleNamesALL, alleles2) #store alleles (including potential stutters
    
    #Obtain number of typed alleles (may include Q-alelle):
    tmp <- rep(0, nAlleles[m])
    typedRefs = unique( c(which(condOrder>0),knownRef,knownRel) ) #get unique of typed referneces
    for(k in typedRefs) { #for each typed refs
      av = unlist(refData[[loc]][[k]])
      if(length(av)==0) next #skip if not found
      av[!av%in%alleles] = Qallele #convert non-observed to Q-allele
        
      ind <- which( alleles%in%av )
      tmp[ind] = tmp[ind] + (length(ind)==1) + 1 #add twice sampled if homozygote 
    } #names(tmp)=alleles
    nTyped <- sum(tmp) #number of total sampled (for each loci)
    maTyped <- tmp  #add vector of typed
    
    #Obtain conditional genotypes; Assign genotypes of known references to knownGind-matrix
    #assign references to knownGind-matrix by values of Glist
    Gmat = numeric() #Obtain genotype index matrix: vectorized upper triangular (1,1),(1,2),...,(1,n),(2,2),(2,3),...,(2,n),....,(n,n)  
    for(i in 1:nAlleles[m]) Gmat = rbind(Gmat, cbind( alleles[rep(i,nAlleles[m] - i + 1)], alleles[i:nAlleles[m]] ))
    
    if(!is.null(condOrder) && any(condOrder>0) && MAXnRefs>0) {
      refDataLoc = refData[[loc]][refNamesAll] #obtain correct order of references
      for(k in seq_along(refNamesAll)) { #for each typed refs
        if(k>length(condOrder) || condOrder[k]==0) next
        av = unlist(refDataLoc[[k]])
        if(length(av)==0) next 
        if(length(av)==1) av = rep(av,2) #impute if exactly one allele
        if(length(av)>2) stop("References can't have more than two alleles!") 
        av[!av%in%alleles] = Qallele #convert non-observed to Q-allele
        Gind1 <- av[1]==Gmat[,1] & av[2]==Gmat[,2]
        Gind2 <- av[2]==Gmat[,1] & av[1]==Gmat[,2]
        GindUse =  which(Gind1 | Gind2) #obtain genottype index to use
        if(length(GindUse)==0) stop(paste0("At marker ",loc,": A reference profile was recorded with a rare allele not in evidence. 
			This error occur when the evidence contains all alleles in the frequency database. Please improve frequency database to enable calculation."))
        knownGind[condOrder[k],m] = GindUse - 1 #subtract with one since we work from 0-indice
      }
    }
    nKnowns[m] = sum(knownGind[,m] > -1)  #obtain number of known
    nUnknowns[m] = nC - nKnowns[m] #obtain number of unknowns for marker
    nGenos[m] = round( nAlleles[m]*(nAlleles[m]+1)/2) #number of genotypes (one contributor)
    nJointGenos[m] = nGenos[m]^nUnknowns[m] #get number of joint combinations
    
    #KINSHIP MODULE
    if(!is.null(knownRel) && ibd0[1]<1 ) {
       av = unlist(refData[[loc]][[knownRel]])
       if(length(av)>0) {
         if(length(av)==1) av = rep(av,2) #impute if exactly one allele
         if(length(av)>2) stop("References can't have more than two alleles!") 
         av[!av%in%alleles] = Qallele #convert non-observed to Q-allele
         Gind1 <- av[1]==Gmat[,1] & av[2]==Gmat[,2]
         Gind2 <- av[2]==Gmat[,1] & av[1]==Gmat[,2]
         relGind[m] = which(Gind1 | Gind2) - 1 #subtract with one since we work from 0-indice
       } 
    }
    ##############################################################
    
    nTypedLong = c(nTypedLong, nTyped) #Obtainined number of typed
    maTypedLong = c(maTypedLong, maTyped) #obtained number of typed
    BWtoLong = c(BWtoLong, BWto) #last index is dummy variable (Dropout doesn't stutter).
    FWtoLong = c(FWtoLong, FWto) #last index is dummy variable (Dropout doesn't stutter). THerefor need QalleleIndex
    BWfromLong = c(BWfromLong, BWfrom) #last index is dummy variable (Dropout doesn't stutter).
    FWfromLong = c(FWfromLong, FWfrom) #last index is dummy variable (Dropout doesn't stutter). THerefor need QalleleIndex
  } #end for each rows
  startIndMarker_nAlleles <- c(0,cumsum(nAlleles)) #get observed alleles start position of marker in long vector
  startIndMarker_nAllelesTot <- c(0,cumsum(nAlleles+nPotStutters)) #get observed alleles start position of marker in long vector
  startIndMarker_nAllelesReps <- c(0,cumsum(nAlleles*nRepMarkers)) #get observed alleles start position of marker in long vector
  startIndMarker_nJointGenos =  c(0,cumsum(nJointGenos)) #get number of genotype outcome start position of marker in long vector
  #Avoid that there are any missing markers for any of the referenecs
  #if(NOK>0 && any(knownGind<0)) stop("Implementation does not support missing reference profiles!")
  if(useDEG) {
    basepairLong = (basepairLong-125)/100 #rescale before
  } else {
    basepairLong = rep(0,sum(nAlleles)) #no base pair information (all zero)
    #length(basepairLong)==sum(nAlleles)
  }

  #Create a warning/note when there is a typed related
  if(any(relGind>=0) && !all(relGind>=0)) warning("WARNING: Missing markers for the typed related should be deselected to avoid wrong LR!") 
  
   #Prepare output
   c = list(nStutterModels=as.integer(nStutterModels),nMarkers=as.integer(nLocs),nAlleles=as.integer(nAlleles),
          startIndMarker_nAlleles=as.integer(startIndMarker_nAlleles),startIndMarker_nAllelesReps=as.integer(startIndMarker_nAllelesReps), 
           peaks=as.numeric(YvecLong),freqs=as.numeric(FvecLong),dropinWeight=as.numeric(DvecLong),
           nTyped=as.numeric(nTypedLong),maTyped=as.numeric(maTypedLong),basepairs=as.numeric(basepairLong),
           BWto=as.integer(BWtoLong),FWto=as.integer(FWtoLong), BWfrom=as.integer(BWfromLong),FWfrom=as.integer(FWfromLong),
           nPotStutters=as.integer(nPotStutters),startIndMarker_nAllelesTot=as.integer(startIndMarker_nAllelesTot),
           nRepMarkers= as.integer(nRepMarkers),nReps = as.integer(nReps), markerNames=locs,QalleleIndex=as.integer(QalleleIndex),
           dropinProb=as.numeric(pCv),fst=as.numeric(fstv), AT=as.numeric(ATv),lambda=as.numeric(lambdav),
           knownGind=as.integer(knownGind),NOK=as.integer(NOK), alleleNames=alleleNames,alleleNamesALL=alleleNamesALL, repNames=repNames,
           ibd=as.numeric(ibd0),relGind=as.integer(relGind), nC=as.integer(nC),nKnowns=as.integer(nKnowns),nUnknowns=as.integer(nUnknowns),
           nJointGenos=as.integer(nJointGenos),startIndMarker_nJointGenos=as.integer(startIndMarker_nJointGenos))
   c$useDEG=useDEG #check of whether to use DEG
   c$genoList=genoList #add genotype list
   c$refNamesCond = refNamesCond #add reference names that are conditoned on (correct order)
   c$hasKinship = ibd0[1]<1 #any(relGind>=0) #check if kinship were defined
   
   #Additional objects required for non-fast version (calcloglik_allcomb/calcloglik_cumprob)
   #Stutters: Only required stutter shifts are stored (in simlar from/to vectors)
   c$nStutters <- rep(0L,nLocs) #init as zero
   c$stuttFromInd <-  c$stuttToInd <- c$stuttParamInd <- integer()
   for(m in 1:nLocs) {
     for(a in 1:nAlleles[m]) { #traverse only "observed alleles"
       aind = startIndMarker_nAlleles[m] + a #obtain allele index
       if(BWS && BWtoLong[aind]>(-1)) { #add BW-stutter contr
         c$stuttFromInd = c(c$stuttFromInd,a-1) #indicate allele index of where STUTTER FROM
         c$stuttToInd = c(c$stuttToInd, BWtoLong[aind]) #indicate allele index of where STUTTER FROM
         c$stuttParamInd = c(c$stuttParamInd,0L) #stutter type "0"
         c$nStutters[m] = c$nStutters[m] + 1L #increment
       }
       if(FWS && FWtoLong[aind]>(-1)) { #add BW-stutter contr
         c$stuttFromInd = c(c$stuttFromInd,a-1) #indicate allele index of where STUTTER FROM
         c$stuttToInd = c(c$stuttToInd,FWtoLong[aind]) #indicate allele index of where STUTTER FROM
         c$stuttParamInd = c(c$stuttParamInd,1L) #stutter type "1"
         c$nStutters[m] = c$nStutters[m] + 1L #increment
       }
     }
   }
   c$stuttFromInd <- as.integer(c$stuttFromInd)
   c$stuttToInd <- as.integer(c$stuttToInd)
   c$startIndMarker_nStutters = as.integer( c(0,cumsum(c$nStutters)) )
  return(c)
}
