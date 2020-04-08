#' @title calcRMPfst
#' @author Oyvind Bleka
#' @description A helpfunction for calculating the random match probability for a person of interest (given fst)
#'
#' @param dat A object retured from getData function
#' @param POIind Index of the references being Person of interest (POI)
#' @param condInd Index of conditional references 
#' @param fst The co-ancestry coefficient (theta-correction). Can be a vector (must contain the marker names)
#' @return ret A vector with random match probabilities for each markers
#' @export

calcRMPfst = function(dat,POIind=1,condInd=NULL,fst=0) {
  if(is.null(dat$refData)) stop("No references found in dat object. Please insert reference data.")
  locs = names(dat$popFreq)
  rmp = rep(1,length(locs)) #a vector for random match probabilities
  names(rmp) = locs
  
  #Check theta correction
  if(length(fst)==1) {
    fstv = rep(fst,length(locs)) #common theta correction for all markers
  } else {
    fstv = setVecRightOrder(fst,  locs)
  }
  
  for(m in 1:length(locs)) { #traverse each loci
    loc = locs[m] #get locus name
    refdat = dat$refData[[loc]] #extract reference data
    numRefs = length(refdat) #get number of refs
    if( POIind>numRefs ) stop(paste0("Marker ",loc,": Specified POI index exceeded number of references."))
    if( !is.null(condInd) && any(condInd>numRefs) ) stop("Some of the conditional indices exceeded the number of references.")
    poi = refdat[[POIind]] #get alleles of poi
    refs = unlist(refdat[condInd]) #get alleles of other references
    tmp = calcGjoint(freq=dat$popFreq[[loc]],nU=1,fst=fstv[m],refK=c(poi,refs)) #get genotype prob conditioning on POI and other refs
    ind = which( (tmp$G[,1]==poi[1] & tmp$G[,2]==poi[2]) | (tmp$G[,1]==poi[2] & tmp$G[,2]==poi[1])) #obtain index to use
    if(length(ind)==0) next #skip if none found
    rmp[m] = tmp$Gprob[ind][1] #insert geno
  }
  return(rmp) #get LR of POI (maximum)
}