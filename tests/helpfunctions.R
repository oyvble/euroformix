
#helpfunction with specified tolerance
expect = function(x,y,tol=testthat_tolerance()) { 
  expect_equal(as.numeric(x),as.numeric(y),tolerance = tol)
}

compareValid = function(x,y,s0=3)  expect(sort(round(x,s0)),sort(round(y,s0)))

#Helpfunction to get likelihood per genotypes (cal in R)
#locs=NULL;modelDEG=TRUE;modelStutt=FALSE
#ibd0=NULL; refRel=NULL
getLogLiki = function(thhat, dat, NOC, cond,pCv,ATv,fstv,lamv, kit0=NULL, ibd0=NULL, refRel=NULL, locs=NULL, modelDEG=TRUE,modelStutt=TRUE) {
  Qallele = "99"
  
  #exctract params:  
  mx=thhat[1:NOC]
  shape0 = 1/(thhat[NOC+2]^2) #get shape param
  scale0 = thhat[NOC+1]/shape0 #get scale param
  
  if(modelDEG) beta =  thhat[NOC+3] #degrad slope param
  if(modelStutt) {
    xiB = thhat[NOC + 4 - as.integer(!modelDEG)] #backwards stutter prop (shift index if no deg model)
    xiF = thhat[NOC + 5 - as.integer(!modelDEG)] #forward stutter prop (shift index if no deg model)
  }
  
  if(is.null(locs)) locs = names(dat$popFreq) #select loci to traverse
  logLikv = rep(NA,length(locs)) #loci to check values for
  names(logLikv) = locs
  for(loc in locs) {
#  loc=locs[6]
    freq = dat$popFreq[[loc]]
    Aset <- Aset0 <- names(freq) #get allele outcome 
    if( Qallele%in%Aset ) Aset0 = head(Aset,-1) #keep only non-dropouts
    
    if(modelStutt) {
      Astutt = c(as.numeric(Aset0)-1,as.numeric(Aset0)+1) #get stutter alleles
      Aset = c(Aset,setdiff(Astutt,Aset0)) #insert potential stutters last
      FWstutt = match(Aset,as.character(as.numeric(Aset)+1)) #alleleinds to stutter from
      BWstutt = match(Aset,as.character(as.numeric(Aset)-1)) #alleleinds to stutter from
    }
    
    #Prepare base pair info for kits:
    if(modelDEG) {
      bpv = rep(NA,length(Aset)) #obtain base pairs for kit
      kitinfo = getKit(kit0) #get kitinfo
      subkit <- subset(kitinfo,toupper(kitinfo$Marker)==loc) #control on dyer
      for(aa in Aset) { #for each allele
        ind <- which(subkit$Allele==aa)
        if(length(ind)==0) ind <- which.min(abs(as.numeric(subkit$Allele) - as.numeric(aa))) #nearest neighbour (in allele name)-> maximum of the alleles
        bpv[aa==Aset] <- subkit$Size[ind] #insert base pair
      }
    }
    
    #Assuming NOK known contribtions
    nRefs = length(dat$refData[[loc]]) #number of refs
    nAG = matrix(0,nrow=length(Aset),ncol=NOC,dimnames = list(Aset)) #create allele counting matrix for unknown contributors
    for(k in seq_len(nRefs)) { #traverse each ref
      if(cond[k]==0) next 
      condRef =  dat$refData[[loc]][[k]] #get alleles of conds
      if(length(condRef)!=2) next #skip if not 2 alleles
      nAG[,cond[k]] = table(factor(condRef,level=Aset)) #insert contribution
    }
    unknownContrs = which(colSums(nAG!=0)==0) #obtain number of unknowns
    nUnknowns = length(unknownContrs)
    
    #Prepare PH data (same order as Aset):
    Ylist = list() #list per samples
    for(sample in names(dat$samples)) {
      Ylist[[sample]] = rep(0,length(Aset))  #preparing PHs after extending allele vector
      Ylist[[sample]][match(dat$samples[[sample]][[loc]]$adata,Aset)] = dat$samples[[sample]][[loc]]$hdata  #insert PHs
    }
    
    nAG0 = nAG #temporary store
    refR = NULL
    if(!is.null(refRel))  refR = unlist( dat$refData[[loc]][refRel] ) #allleles of releated
#    nU=max(1,nUnknowns);fst=fstv[loc];refK=unlist( dat$refData[[loc]] );ibd = ibd0;sortComb = FALSE
    Glist = calcGjoint(freq=freq,nU=max(1,nUnknowns),fst=fstv[loc],refK=unlist( dat$refData[[loc]] ),refR=refR,ibd = ibd0,sortComb = FALSE) #include all typed refs here
    Gset = Glist$G #get allele out come of unknowns
    
    pEvid = 0 #obtain probabiliity of evidence (sum over all genotypes)
    for(gind2 in 1:nrow(Gset)) { #traverse all genorypes (C2)
      for(gind1 in 1:nrow(Gset)) { #traverse all genorypes (C1)
#        gind1 <- gind2 <- 1
        nAG = nAG0 #copy back
        if(nUnknowns==1) {
          nAG[,unknownContrs] = table(factor(Gset[gind1,],level=Aset)) #insert unknown contribution
        } else if(nUnknowns==2) {
          nAG[,unknownContrs[1]] = table(factor(Gset[gind1,],level=Aset)) #insert unknown contribution
          nAG[,unknownContrs[2]] = table(factor(Gset[gind2,],level=Aset)) #insert unknown contribution
        } else if(nUnknowns>2) {
          stop("Script not supporting more than 2 unknown")
        } 
        mui <-  c(nAG[,drop=FALSE]%*%mx)*shape0 #expected contribution
        if(modelDEG) mui  <- mui*exp( (bpv-125)/100*log(beta) )  #in case of degradation
        
        #Delegate stutterprod
        if(modelStutt) {
          stuttB <- mui[BWstutt]*xiB #backward-stutter parts
          stuttF <- mui[FWstutt]*xiF #forward-stutter parts
          indBW = !is.na(stuttB)
          indFW = !is.na(stuttF)
          
          indLooseStutt = indBW & indFW #index of alleles which are assumed to not loose stutter product
          mui[indLooseStutt] = mui[indLooseStutt]*(1- (xiB+xiF)) #loose stutter products
          mui[indBW] = mui[indBW] + stuttB[indBW]
          mui[indFW] = mui[indFW] + stuttF[indFW]
        }
        
        vali = 0 
        for(sample in names(dat$samples)) { #traversing all samples
          #Divide set into dropin, contr and dropout
          psiDI <- which( Ylist[[sample]]>0 & mui==0 ) #drop-in (not a stutter). ANY DROPIN?
          psiYmu <- which( Ylist[[sample]]>0 & mui>0 )  #contributing to model and observed PH
          psiDO <- which( Ylist[[sample]]==0 & mui>0 )  #dropout elem
          
          if(length(psiYmu)>0) vali =  vali + sum(dgamma(Ylist[[sample]][psiYmu],shape=mui[psiYmu],scale=scale0,log=TRUE))   # CONTRIBUTION TO PH
          if(length(psiDO)>0) vali = vali + sum(pgamma(ATv[loc],shape=mui[psiDO],scale=scale0,log=TRUE)) #CONTR TO DROPOUT
          if(length(psiDI)>0) vali = vali + sum( dexp( Ylist[[sample]][psiDI]-ATv[loc], rate=lamv[loc], log=TRUE) + log(pCv[loc]*freq[psiDI]) ) #CONTR TO DROPIN: Scaled for each dropin. THIS IS OK
          if(length(psiDI)==0) vali = vali + log(1-pCv[loc]) #in case of no dropin
        }
        genProb = 1 #INSERT Prob genotypes
        if(nUnknowns==1) genProb = Glist$Gprob[gind1]
        if(nUnknowns==2) genProb = Glist$Gprob[gind1,gind2] #related is last
        
        pEvid = pEvid + exp( vali + log(genProb)) #sum up
        if(nUnknowns==0) break #stop loop if no unknowns
      } #end for each genotype (C1)
      if(nUnknowns < 2) break #stop loop if less than 2 unknown
    } #end for each genotype (C2)
    logLikv[loc] = log(pEvid)
  } #end for each locus
  return(logLikv)
} 

#Helpfunction for checking qualitative model after model fit (iterating over 1 genotype set)
checkQualModel = function(mle) {
  
  #COMPARE PER MARKER RESULTS WITH MANUAL DERIVED:
  logLikv = mle$logliki
  expect_equal(sum(logLikv),mle$loglik) #checking joint loglik
  
  thhat = mle$pDhatContr #obtain drop-out probs per contributor
  cond = mle$model$condOrder
  NOC = mle$model$nC
  pCv = mle$model$prC
  fstv = mle$model$fst
  
  #Obtain data  
  popFreq = mle$model$popFreq
  refData = mle$model$refData
  samples = mle$model$samples
  
  #exctract params:  
  locsCheck = names(logLikv) #loci to check values for
  locsVal = setNames(rep(NA,length(locsCheck)),locsCheck)
  for(loc in locsCheck) {
    #loc=locsCheck[1]
    freq = popFreq[[loc]]
    Aset = names(freq) #get allele outcome 
    
    #Assuming a given number of known contribtions
    condUseIdx = which(cond>0)
    nAG = matrix(0,nrow=length(Aset),ncol=NOC,dimnames = list(Aset)) #create allele counting matrix for unknown contributors
    for(ind in condUseIdx) nAG[,cond[ind]] = table(factor(refData[[loc]][[ind]],level=Aset)) #insert contribution
    
    #Obtain all combinations for unknown
    unknownIdx = which(colSums(nAG>0)==0) #get index of unknowns
    nU = length(unknownIdx)
    nG = 1
    if(nU>0) {
      Glist = calcGjoint(freq=freq,nU=nU,fst=fstv[loc],refK=unlist(refData[[loc]])) #don't use [which(cond>0)] since we include all
      Gset = Glist$G #get allele out come of unknowns
      gProb = log(Glist$Gprob)
      nG1 = nrow(Gset) #number of genotypes (1p) 
      nG = length(gProb) #number of joint genotypes (Up)
    }
    #Traversing all combinations for unknown
    pEvid = 0 #obtain probabiliity of evidence (sum over all genotypes)
    pGeno = 0
    nAG0 = nAG #temporary store
    for(gind in seq_len(nG) ) { #traverse all genorypes
      nAG = nAG0
      if(nU>0) {
        GsetInd = rep(0,nU)
        genoPower = 1
        for(u in unknownIdx) {
          #u = 1
          GsetInd[which(unknownIdx==u)] <- g1ind <- as.integer( ((gind-1)/genoPower) %% nG1) + 1
          genoPower = genoPower*nG1
          nAG[,u] = nAG[,u] + table(factor(Gset[g1ind,],level=Aset)) #insert contribution
        }
        pGeno = gProb[rbind(GsetInd)]
      }
      
      #temporary caluclation
      pDprodTmp = 1 #calculte temporary
      for(k in 1:NOC) { #traverse each contributor
        pDprodTmp = pDprodTmp*thhat[k]^nAG[,k]
      } 
      mui = 1-pDprodTmp #rowSums(nAG) #obtain per-allele contribution
      
      #calculate likelihood for observations
      vali = 0 
      for(sample in names(samples)) { #for each replicates
        #Divide set into dropin, contr and dropout
        #      sample=names(samples)[1]
        yvec = rep(0,length(Aset))  #preparing PHs after extending allele vector
        yvec[match(samples[[sample]][[loc]]$adata,Aset)] = samples[[sample]][[loc]]$hdata  #insert PHs
        
        psiDI <- which( yvec>0 & mui==0 ) #drop-in (not a stutter). ANY DROPIN?
        psiYmu <- which( yvec>0 & mui>0 )  #contributing to model and observed PH
        psiDO <- which( yvec==0 & mui>0 )  #dropout elem
        
        if(length(psiYmu)>0) vali =  vali + sum(log(1-pDprodTmp[psiYmu]))   # CONTRIBUTION TO PH
        if(length(psiDO)>0) vali = vali + sum(log(pDprodTmp[psiDO])) #CONTR TO DROPOUT
        if(length(psiDI)>0) {
          vali = vali + sum(log(pCv[loc]*freq[psiDI])) #CONTR TO DROPIN: Scaled for each dropin. THIS IS OK
        } else {
          vali = vali + log(1-pCv[loc]) #in case of no dropin
        }
      } #end for each sample
      pEvid = pEvid + exp( vali + pGeno) #sum up
    } #end for each genotype
    locsVal[loc] = log(pEvid) #return()
    expect_equal(as.numeric(logLikv[loc]),as.numeric(log(pEvid)))
  }
  return(locsVal)
}



#helpfunction to obtain cumulative probs using integration
#THIS VERSION DOES NOT (YET) WORK
getValidProbs = function(mle) {
  library(cubature)
  
  maxYobs <- max(sapply(mle$model$samples,function(x) sapply(x,function(y) max(y$hdata)))) #max observation
  nAtotObs = sum(sapply(mle$model$samples,function(x) sapply(x,function(y) length(y$adata)))) #number of observed alllees
  thhat = mle$fit$thetahat2[c$nC+1:2]
  
  #Obtain large Y:
  alphaQ <- 0.001 #ensure very far out in quantile (used for estimating probs in gamma-distribution).
  alpha2 <- alphaQ/nAtotObs #"bonferroni outlier"
  maxYexp <- qgamma(1-alpha2,2/thhat[2]^2,scale=thhat[1]*thhat[2]^2) #max observation in theory
  maxY <- ceiling(max(maxYobs,maxYexp)) #get max observed
  
  
  #Use logLiki and integrals to obtain 
  mle2 = mle #copy
  c2 <- c <- mle2$prepareC
  locs = c$locs
  #Following code taken from logLiki
  loglikv <- loglikv2 <- logLiki(mle)
  
  #Step 1) Calculate L(E|thetahat) for each marker
  startIndPS = 0; #start marker index for potential stutters (own vectors)
  startIndMarker1=0;#start marker index (1 rep)
  startIndMarker2=0;#start marker index (nRep[m] reps)
  cumProb = NULL
  for(m in seq_along(locs)) {
    loc = locs[m]
    ind1 = startIndMarker1 + 1:c$nA[m] #get index of 1 repitition
    ind2 = startIndMarker2 + 1:(c$nA[m]*c$nRep[m]) #get index of 1 repitition
    indPS = startIndPS +  1:c$nPS[m] #get index of potential stutters
    if(c$nPS[m]==0) indPS = numeric()
    
    genoinds = (c$nC*(m-1)+1):(c$nC*m) #get index of known genotype index (knownGind and relGind)
    
    c2$nM = 1 #restrict to only 1 marker
    c2$nReps = c$nReps[m]
    c2$nA = c$nA[m]
    c2$nPS <- c$nPS[m]
    c2$nTypedLong = c$nTypedLong[m]
    c2$NOK = c$NOK[m]
    c2$FvecLong = c$FvecLong[ind1]
    c2$maTypedLong = c$maTypedLong[ind1]
    c2$basepairLong = c$basepairLong[ind1]
    c2$BWvecLong = c$BWvecLong[ind1]
    c2$FWvecLong = c$FWvecLong[ind1]
    c2$BWPvecLong = c$BWPvecLong[indPS]
    c2$FWPvecLong = c$FWPvecLong[indPS]
    c2$knownGind = c$knownGind[genoinds] 
    c2$relGind = c$relGind[genoinds]
    c2$locs =  c$locs[m]
    
    mle2$model$popFreq = mle$model$popFreq[loc]
    mle2$model$prC = mle$model$prC[[loc]]
    mle2$model$fst = mle$model$fst[[loc]]
    mle2$model$lambda = mle$model$lambda[[loc]]
    mle2$model$threshT = mle$model$threshT[[loc]]

    yobs = c$YvecLong[ind2] #obtain observed vals
    
    #CHECK:
    c2$YvecLong = yobs
    mle2$prepareC = c2 #insert modified object
    loglikv2[m] = logLiki(mle2) 

    #Perform integration for each allelel
    for(a in which(yobs>0)) { #traverse only where PH observed
      fun_Ya = function(x) {
        c2$YvecLong[a] = x #insert PH
        mle2$prepareC = c2 #update object
        exp(logLiki(mle2)) #take exponent
      }

      AT = mle2$model$threshT #get threshold
      minY = AT-1 #lower limit
      #1) Integrate from max 
      tol=0.0000001
      num = cubature::adaptIntegrate(Vectorize(fun_Ya),minY,yobs[a],tol=tol)[[1]]
      denom = cubature::adaptIntegrate(Vectorize(fun_Ya),minY,maxY,tol=tol)[[1]]
      cumProb = c(cumProb, num/denom)
    }
        
    #Update indices for next marker:
    startIndMarker1 = startIndMarker1 + c$nA[m] #get start position of marker m+1 (1 rep)
    startIndMarker2 = startIndMarker2 + c$nA[m]*c$nRep[m]; #get start position of marker m+1 (nRep), used only for Peaks only
    startIndPS = startIndPS + c$nPS[m]; #add number of potential stutters
  }
  return(cumProb)
}





