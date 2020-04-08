#' @title prepareC 
#' @author Oyvind Bleka
#' @description prepareC is used in the functions contLikMLE,contLikMarg and contLikMCMC to prepare input to C-call
#' @details Assumes that all PH below threshold are removed. The function builds the data input to the C-code
#' @param nC Number of contributors in model.
#' @param samples A List with samples which for each samples has locus-list elements with list elements adata and hdata. 'adata' is a qualitative (allele) data vector and 'hdata' is a quantitative (peak heights) data vector.
#' @param popFreq A list of allele frequencies for a given population (possibly including "99" Q-alleles).
#' @param refData Reference objects with list element [[s]][[i]]. The list element has reference-list with list-element 's' having a loci-list adata with list-element 'i storing qualitative data.
#' @param condOrder Specify conditioning references from refData (must be consistent order). For instance condOrder=(0,2,1,0) means that we restrict the model such that Ref2 and Ref3 are respectively conditioned as 2. contributor and 1. contributor in the model. 
#' @param knownRef Specify known non-contributing references from refData (indices). For instance knownRef=(1,2) means that reference 1 and 2 is known non-contributor in the hypothesis. This affectes coancestry correction.
#' @param kit shortname of kit: Obtained from getKit()
#' @param knownRel Specify the index of the related contributing reference from refData (one index). For instance knownRel=2 means that unknown1 is related to reference 2 with ibd specified relationship.
#' @param ibd the identical by decent coefficients list of the relationship for each unknowns (specifies the type of relationship). Default is NULL, meaning no related inds
#' @param fst The co-ancestry coefficient. Default is 0. Can be a vector following markers.
#' @param incS A boolean whether potential BW stutters are included
#' @param incFS A boolean whether potential FW stutters are included
#' @return ret A list of data input to call the C-code with
#' @export 
#' @examples
#' \dontrun{
#' kit = "ESX17"
#' popfn = paste(path.package("euroformix"),"tutorialdata","FreqDatabases",
#'  paste0(kit,"_Norway.csv"),sep=.Platform$file.sep)
#' evidfn = paste(path.package("euroformix"),"examples",paste0(kit,"_3p.csv"),
#'  sep=.Platform$file.sep)
#' reffn = paste(path.package("euroformix"),"examples",paste0(kit,"_refs.csv"),
#'  sep=.Platform$file.sep)
#' popFreq = freqImport(popfn)[[1]] #obtain list with population frequencies
#' samples = sample_tableToList(tableReader(evidfn))
#' dat = prepareData(samples,popFreq=popFreq) #obtain data to use for analysis
#' prepCobj = prepareC(nC=2,samples=dat$samples,popFreq= dat$popFreq,kit=kit)
#' }


prepareC = function(nC,samples,popFreq,refData=NULL,condOrder=NULL,knownRef=NULL,kit=NULL,knownRel=NULL,ibd=NULL,fst=0,incS=FALSE,incFS=FALSE){
 Qallele="99" #Name of allele given if missing in evidence. Defualt is 99. This is important when considering the degradation model since 99 is closest to maximum allelein a locus. 
 LUSsymbol="_" #a character symbol used to separate repeatunit and LUS.

 #CHECK INPUT OF Relatedness elements:
 if(!is.null(knownRel) && length(knownRel)>1 ) stop("Not implemented: Multiple related not possible!")
 if(!is.null(ibd) && length(ibd)!=3 ) stop("Not implemented: Multiple ibs elements not possible! Ibs must be a vector")

  #Note: Supports Invariant order of markers and caseletters!!
 #Supports replicates in the samples. 
 #Evaluates markers which are in both mixture and frequency (hence number of markers)
 nS <- length(samples) #number of replicates (assume same number for each reps)
 #test 1) Require same markers in all replicates. 
 locs <- lapply(samples,function(x) toupper(names(x)))  #marker names for each replicates
 nLs <- sapply(locs,length) #number of markers per replicaeetes
 nL <- unique(nLs) #get unique number of markers
 if(length(nL)>1) stop("Number of markers in the replicates was not the same")
 locs <- unique(unlist(locs)) #get unique markernames
 if(nL!=length(locs)) stop("Different markers was specified in the replicates")
 names(popFreq) <- toupper(names(popFreq)) #toupper case!
 locs <- intersect(locs,names(popFreq)) #take intercept with locus in evidence and population-freq
 #print( paste0("Evaluated loci: ", paste0(locs,collapse=",") ) )

 #convertion of values in popFreq, mixData and Glist$G:
 #loci-order follows as in mixData: "locs". Rearrange names:
 popFreq <- popFreq[locs] #order popFreq to mixData-order
 
 #Get list of genotypes for each markers (MUST BE SAME ORDER AS THOSE CREATED IN C++ code!!)
 Gset <- list() 
 for(loc in locs) Gset[[loc]] <- calcGjoint(freq=popFreq[[loc]])$G #get genotypes


 #Fix references as known contributors: Assign genotypes of known references to knownGind-matrix
 NOK = rep(0,nL) #number of known contributors per loci
 knownGind <- matrix(-1,ncol=nL,nrow=nC) #default is no references (=-1)
 #assign references to knownGind-matrix by values of Glist
 if(!is.null(refData) && !is.null(condOrder) && any(condOrder>0)) {
   names(refData) <- toupper(names(refData)) #set marker names to toupper case!
   for(loc in locs) {
     locind = which(loc==locs)
     subRef <- refData[[loc]] #take out relevant reference
     if(length(subRef)==0)  stop(paste('Missing locus (',loc,') in refData.',sep=''))
     for(k in 1:length(subRef)) { #for each reference
       if(length(subRef[[k]])==0) next #updated from v.0.6.4: Allowing unknown contributors for missing markers.
       if(length(subRef[[k]])!=2) stop("References need to have exactly two alleles.") #updated from v.0.6.4
       if(condOrder[k]>0) {
         NOK[locind] = NOK[locind] + 1 #add known contributor
         Gind1 <- subRef[[k]][1]==Gset[[loc]][,1] & subRef[[k]][2]==Gset[[loc]][,2]
         Gind2 <- subRef[[k]][2]==Gset[[loc]][,1] & subRef[[k]][1]==Gset[[loc]][,2]
         knownGind[condOrder[k],locind] = which(Gind1 | Gind2) - 1 #subtract with one since we work from 0-indice
       }
     }
   }
 }
 #NOU = nC - NOK #number of unknowns per markers
 startUindex = max(NOK) + 1 #index of first unknown to be related
 
 #Assign genotypes of related references to relGind-matrix
 relGind <- matrix(-1,ncol=nL,nrow=nC) #default is no references (=-1)
 #assign references to knownGind-matrix by values of Glist
 if(!is.null(refData) && !is.null(knownRel) ) {
   names(refData) <- toupper(names(refData)) #set marker names to toupper case!
   for(loc in locs) {
     locind = which(loc==locs)
     subRef <- refData[[loc]] #take out relevant reference
     if(length(subRef)==0)  stop(paste0("Missing locus ",loc," in refData."))
     for(k in knownRel) { #for each related reference
       if(length(subRef[[k]])==0) {
         print(paste0("Missing alleles at locus ",loc," for related reference ",names(subRef)[k]))
         next
       }
       if(length(subRef[[k]])!=2) stop("References need to have exactly two alleles.") #updated from v.0.6.4
       
       Gind1 <- subRef[[k]][1]==Gset[[loc]][,1] & subRef[[k]][2]==Gset[[loc]][,2]
       Gind2 <- subRef[[k]][2]==Gset[[loc]][,1] & subRef[[k]][1]==Gset[[loc]][,2]
       relGind[startUindex,locind] = which(Gind1 | Gind2) - 1 #subtract with one since we work from 0-indice
     }
   }
 }
 anyRel =  any(relGind!=-1) #any related references (at least one marker)?
 ibd0 =  c(1,0,0) #unrelatedness is default
 ibdList = list()
 for(cc in 1:nC) ibdList[[cc]] =   ibd0 #insert default values for each contributors  
 if(!is.null(ibd) && anyRel && ibd[1]<1) { #critetion for relatedness
   ibdList[[startUindex]] = ibd #insert given ibs to the unknown
 } else {
   anyRel = FALSE #no related otherwise
 }
 ibdLong = unlist(ibdList) #long vector over all contributors

 #INITS:
 nM = length(locs) #[1] #num markers to evaluate with
 nA <- nPS <- rep(0,nM) #num alleles and extra potential stutters  
 nReps <- rep(nS,nM) #num replicates per markers (assume same as nS for now)
 YvecLong <- FvecLong <- numeric()
 BWvecLong <- FWvecLong<- numeric() #backward/forward stutter info
 BWPvecLong <- FWPvecLong<- numeric() #potential backward/forward stutter info (in addition)
 basepairLong <- numeric() #base pair information
 if(!is.null(kit)) {
   kitinfo = euroformix::getKit(kit) #get kit information
 } 
 avL = list() #store used alleles per marker
 
 #HELPFUNCTION USED FOR LUS:
 getLUS_BWstutter = function(x) { #get the BW stuttered variants of LUS
   tmp = as.numeric(x) 
   stutt2 <- paste0(tmp[1]-1,LUSsymbol,tmp[2]-1) #BW stutter-allele
   if(length(tmp)>2) stutt2 <- paste0(c(stutt2,tmp[-(1:2)]),collapse=LUSsymbol) #add other LUS variants
   return(stutt2)
 }

 getLUS_FWstutter = function(x) { #get the FW stuttered variants of LUS
   tmp = as.numeric(x) 
   stutt2 <- paste0(tmp[1]+1,LUSsymbol,tmp[2]+1) #FW stutter-allele
   if(length(tmp)>2) stutt2 <- paste0(c(stutt2,tmp[-(1:2)]),collapse=LUSsymbol) #add other LUS variants
   return(stutt2)
 }
 
 #TRAVERSE FOR EACH MARKER:
 #PREPARE FREQUENCY AND PH INFO
 for(m in 1:nM) { #m=7
   loc = locs[m] #for selected loc (already upper)
   freq = popFreq[[loc]] #get allele freqs
   avL[[loc]] = names(freq) #get orignal allele outcome (including Qallele) 
   av = setdiff(names(freq),Qallele) #assumed order as in freqs (assume alleles already added!)
   FvecLong = c(FvecLong , freq) #add freqs (Q-allele already included)
      
   #Check if alleles are in LUS format (i.e. contains the _ separation symbol):   
   isLUS <- FALSE
   if(length(av)>0) {
     tmp = grepl(LUSsymbol,av)
     if(all(tmp)) { #if all alleles was with LUS syntax
       isLUS = TRUE  #LUS IS CONSIDERED FOR THE MARKER!
     } else if(any(tmp)) { #if only some alleles are in LUS
       stop(paste0("In marker ",loc,": Only some alleles contained the LUS syntax. All or no alleles must contain this syntax for a marker. Program stops!"))
     }
   }

   #Fix PH vector:      
   nA[m] <- length(freq) #number of allees to travers thgotu 
   locDat = lapply(samples,function(x) x[which(toupper(names(x))==loc)][[1]]) #obtain replicate info
   yv = matrix(0,ncol=nA[m],nrow=nS) #create PH-matrix
   for(r in 1:nS) { #for each replicates (following handle unordered loci under each sample)
     yv[r, match(locDat[[r]]$adata,av)] = locDat[[r]]$hdata #insert PHs
   }
   YvecLong = c(YvecLong , unlist(yv)) #add PHs for obs alleles
   
   if(!is.null(kit)) { #Extracting basepairs
     sub = kitinfo[toupper(kitinfo$Marker)==loc,,drop=FALSE] 
     av0 = av #copy alleles
     if(isLUS) av0 <- as.numeric( sapply(strsplit(av0 ,LUSsymbol),function(x) x[1] ) ) #extracting first allele 
     
     bp = sub$Size[match(as.character(av0),sub$Allele)] #corresponding bp of alleles
     if(any(is.na(bp))) { #if allele not found
       for(j in which(is.na(bp)) ) bp[j] =  sub$Size[which.min(abs(as.numeric(sub$Allele) - as.numeric(av0[j])))] #the closest allele is used
     }
     basepairLong = c(basepairLong,bp) 
     if(Qallele%in%names(freq)) basepairLong = c(basepairLong,max(sub$Size)) #set max size for Q-allele if considered
   } 
   
   #Prepare BWstutter/FWstutter relations:
   suppressWarnings({ av2 = as.numeric(av)}) #convert to numeric
   isNum = ifelse(!any(is.na(av2)),TRUE,FALSE ) #check if all alleles are numeric
   #isNum = all(is.numeric(av)) #check if all alleles are num (stutters may be considered)
   
   nPS[m] = 0  #number of potential stutters (Alleles not part of av)
   if(isNum || isLUS ) {
     if(isNum) {
       BWstutt <- as.character(as.numeric(av)-1) #get corresponding BW stutters
       FWstutt <- as.character(as.numeric(av)+1) #get corresponding FW stutters
     } else { #if it was LUS
       avLUSlist = strsplit(av ,LUSsymbol) #split up elements with LUS separator
       BWstutt <- sapply(avLUSlist, getLUS_BWstutter) #get corresponding BW stutters
       FWstutt <- sapply(avLUSlist, getLUS_FWstutter) #get corresponding FW stutters
     }
     BWind <- match(as.character(av),as.character(BWstutt),nomatch=0)-1 #get index of what allele it receive stutter from  (Allele 1=index0)
     FWind <- match(as.character(av),as.character(FWstutt),nomatch=0)-1 #get index of what allele it receive stutter from  (Allele 1=index0)
     
     Pstutt = numeric() #potential stutters (alleles not in genotype set)
     if(incS) Pstutt = c(Pstutt,BWstutt) #add potential stutters for BW
     if(incFS) Pstutt = c(Pstutt,FWstutt) #add potential stutters for FW
     Pstutt = setdiff(unique(Pstutt),as.character(av) ) #get all potntial non-observed stutters (both FW/BW)
     BWPind <- match(as.character(Pstutt),as.character(BWstutt),nomatch=0)-1 #get corresponding potential BW stutters
     FWPind <- match(as.character(Pstutt),as.character(FWstutt),nomatch=0)-1 #get corresponding potential FW stutters
     nPS[m] <- length(Pstutt) #number of potential stutters
   } else { #stutters not possible to extract
     BWind <- FWind <- rep(-1,nA[m]) #add emtpy array
     BWPind <- FWPind <- numeric() #empty vector
   }
   BWvecLong = c(BWvecLong, BWind, -1) #last index is dummy variable (Dropout doesn't stutter)
   FWvecLong = c(FWvecLong, FWind, -1) #last index is dummy variable (Dropout doesn't stutter)
   BWPvecLong = c(BWPvecLong, BWPind) #add stuttered index of potential stutters
   FWPvecLong = c(FWPvecLong, FWPind) #add stuttered index of potential stutters
 } #end for each rows

 if(is.null(kit)) {
   basepairLong = rep(0,sum(nA)) #no base pair information (all zero)
 } else {
   length(basepairLong)==sum(nA)
   basepairLong = (basepairLong-125)/100 #rescale before
 }
 
 #Get number of typed allelse (must count alleles)
 nTypedLong <- maTypedLong <- integer()
 for(loc in locs) { #for each loci
   av = avL[[loc]] #extract alleles
   tmp <- rep(0, length(av))
   if(!is.null(refData)) {
     typedRefs = unique( c(which(condOrder>0),knownRef,knownRel) ) #get unique referneces
     for(k in typedRefs) { #for each typed refs
       ind <- which( av%in%refData[[loc]][[k]] )
       tmp[ind] = tmp[ind] + (length(ind)==1) + 1 #add twice sampled if homozygote 
     }
   }
   nTypedLong <- c(nTypedLong,sum(tmp)) #number of total sampled (for each loci)
   maTypedLong <- c(maTypedLong,tmp)  #add vector of typed
 }
# names(maTypedLong) = unlist(avL)
 
 retlist = list(nC=as.integer(nC),nReps=as.integer(nReps), nM=as.integer(nM),nA=as.integer(nA),YvecLong=as.numeric(YvecLong),FvecLong=as.numeric(FvecLong),nTypedLong=as.numeric(nTypedLong),maTypedLong=as.numeric(maTypedLong),basepairLong=as.numeric(basepairLong),BWvecLong=as.integer(BWvecLong),FWvecLong=as.integer(FWvecLong),nPS=as.integer(nPS),BWPvecLong=as.integer(BWPvecLong),FWPvecLong=as.integer(FWPvecLong),knownGind=as.integer(knownGind),NOK=as.integer(NOK),locs=locs, anyRel=as.integer(anyRel),relGind=as.integer(relGind), ibdLong=as.numeric(ibdLong) )
 
 return(retlist)
} #end function

