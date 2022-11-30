#' @title calcLRupper
#' @description Obtaining upper boundary LR for POI 
#' @details The function takes a fitted MLE object under Hd (where hypothesis$knownRef is the POI)
#' 
#' @param POIidx Index of POI (in refData)
#' @param mle Fitted object using calcMLE (under Hd)
#' @param scale Whether to scale the RMP when having more than 2 unknowns (recommended)
#' @return Returning theoretical upper boundary LR (log10 scale)
#' @export

calcLRupper = function(POIidx,mle, scale=TRUE) {
  #print(POIidx)
  #mle <<- mle
  model = mle$model
  fst0 = model$fst
  if(length(POIidx)!=1) return(NA) #print("Couldn't calculate upper boundary LR. Returning...")
  
  c = mle$prepareC #obtain stored C-object
  locs = c$markerNames #get evaluating markers

  #$knownGind
  nKnownRefs = c$NOK #number of known refs
  #indKnownGenoCidx = which(model$condOrder==POIidx) #handle if multiple known references
  
  rmp = setNames( rep(1,length(locs)) , locs) #calculated Random match probability
  for(m in seq_along(locs)) {
  #m=1  
    loc = locs[m]
    SImarker = c$startIndMarker_nAlleles[m] #start index of marker
    alleleRange = SImarker + seq_len(c$nAlleles[m]) #obtain allele range
    alleles = c$alleleNames[alleleRange] #obtain allele names
    geno = c$genoList[[loc]] #genotype outcome
    genoIdx = c$knownGind[ nKnownRefs*(m-1) + POIidx ] + 1 #get genotype index (REMEMBER TO ADD +1 to have R index)
    if(genoIdx == 0) next #Skip if missing reference (don't scale on this allele!)
    alleleIdx = match(geno[genoIdx,],alleles) #obtain allele index of POI
    freqs = c$freqs[alleleRange] #allele frequencies
    maTyped = c$maTyped[alleleRange] #number of typed alleles
    totTyped = c$nTyped[m] #number of typed alleles
    fst0 = c$fst[m]

    gprob = 1 #prod(freqRef) #genotype prob
    for(a in 1:2) {
      aind = alleleIdx[a] #get correct index
      gprob = gprob * (fst0*maTyped[aind] + (1-fst0)*freqs[aind]) / (1 + (totTyped-1)*fst0); 
      maTyped[aind] = maTyped[aind] + 1; #update allele count 
      totTyped = totTyped + 1; #update total count
    }
    if(alleleIdx[1]!=alleleIdx[2]) gprob = 2*gprob #get het situation
    rmp[m] = gprob
  }
  
  #THIS LAST BLOCK IS TO TAKE INTO ACCOUNT THE POSSIBILITY THAT 
  #THE LR CAN EXCEED 1/RMP WHEN HAVING MORE THAN 2 unknowns
  scale0 = 1 #default is no scaling
  fst = c$fst
  nCond = c$nKnowns - 1 #Number of conditionals under Hd (handles missing markers here)
  nU = model$nU + 1  #Number of unknowns under Hd (adjust by 1)
  if(nU>=2 && scale) { #if at least 2 unknowns (takes into account that profile corresponding to POI could be a clear Major)
    scale0 = (1+(3+2*nCond)*fst)/(1+(1+2*nCond)*fst)*(1+(4+2*nCond)*fst)/(1+(2+2*nCond)*fst) 
    scale0[rmp==1] = 1 #ensure that scaling is 1 for non-informative markers (CORRECT?)
  }
  upperLR = -sum(log10(rmp/scale0))  #get upper limit of LR (log10 scale), fst taken into account
  return(upperLR)
}

