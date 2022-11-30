#ALGORITHM TO CALCULATE THE RESTRICTED LIKELIHOOD FAST (SYMMETRY IN CONTriBUTORS)

likQualR = function(evids,conds, knowns, nUnknowns, fst, pD, pC, freq) {
 #evids,conds,knowns,nUnknowns,fst, rep(pD,nContr),rep(pD,nContr)^2,pC,freq)

  alleleNames = names(freq) #obtain allele names
  m_nAlleles = length(freq) #number of alles
  m_nGenos = m_nAlleles*(m_nAlleles+1)/2
  
  m_contrMat = matrix(0,nrow=m_nGenos,ncol=m_nAlleles)
  m_alleleMat = matrix(0,nrow=m_nGenos,ncol=2)
  cc = 1
  for(i in seq_len(m_nAlleles)) {
    for(j in i:m_nAlleles) {
      m_contrMat[cc,i] =  m_contrMat[cc,i] + 1
      m_contrMat[cc,j] =  m_contrMat[cc,j] + 1
      m_alleleMat[cc,1] = i
      m_alleleMat[cc,2] = j
      cc = cc + 1
    }
  }
  
  #CREATE A BOOLEAN OF WHETHER EVIDENCE HAS ALLELE IN FREQUENCY DATA
  hasEvidMat = NULL
  hasEvidVec0 = rep(0L,m_nAlleles)
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
  nReplicates = ncol(hasEvidMat) #number of replictes
  
  #Prepared typed info:
  m_nTyped <- m_contrAlleles <- rep(0,m_nAlleles) #number of typed alleles (Each type)
  for(i in seq_len(m_nAlleles)) {
    if(!is.null(conds)) {
      numAlleles = sum(conds==alleleNames[i])
      if(numAlleles>0) {
        m_nTyped[i] = m_nTyped[i] + numAlleles
        m_contrAlleles[i] = m_contrAlleles[i] + numAlleles
      }
    }
    if(!is.null(knowns)) m_nTyped[i] = m_nTyped[i] + sum(knowns==alleleNames[i])
  }
  
  #START RECURSION:
  #genoJointIdx vector with contribution
  m_genoCombConstant = lgamma(nUnknowns+1) #calculate genotype combination (numerator)

  recfun = function(contrIdx, genoJointIdx_rec, genoProb_rec, contrAlleles_rec, nTyped_rec) {

    #traversing all genotypes
    for(genoIdx in seq_len(m_nGenos) ) { 
      
      #RESTRICTION: ORDERED GENOTYPES
      #if( contrIdx>1 && any( genoIdx > genoJointIdx_rec[1:(contrIdx-1)]) ) next 
      if( contrIdx>1 && genoIdx > genoJointIdx_rec[contrIdx-1]) break; 
      
      #Re-copy variables for each new genotype
      genoJointIdx = genoJointIdx_rec
      genoProb = genoProb_rec
      contrAlleles = contrAlleles_rec
      nTyped = nTyped_rec
      
      #Update variables
      genoJointIdx[contrIdx] = genoIdx #insert genotype contribution
      contrAlleles = contrAlleles + m_contrMat[genoIdx,] #add contribtion
      
      for(i in 1:2) { #traverse both alleles in genotype
        alleleIdx = m_alleleMat[genoIdx,i] #obtain allele idx for genotype
        genoProb =  genoProb*(nTyped[alleleIdx]*fst + (1-fst)*freq[alleleIdx]) / (1+(sum(nTyped)-1)*fst)
        nTyped[alleleIdx] = nTyped[alleleIdx] + 1 #add
      }
      if(m_alleleMat[genoIdx,1]!=m_alleleMat[genoIdx,2]) genoProb = 2*genoProb

      #Calculate likelihood if last contribution
      if(contrIdx==nUnknowns) { 
        hasContr = contrAlleles > 0 #bool of which has contributin
        hasNoContr = !hasContr
        dropAlleleContr = pD^contrAlleles#[hasContr] #calculate log  pD^n_a (used two places)
          
        likEvid = 1
        for(r in seq_len(nReplicates)) {
          hasEvid = hasEvidMat[,r]
          
          A_contrSet = which(hasEvid & hasContr) #Contribution set
          B_dropoutSet = which(!hasEvid & hasContr) #Dropout set
          C_dropinSet = which(hasEvid & hasNoContr) #Dropin set
          if(length(A_contrSet)>0) likEvid = likEvid*prod(1-dropAlleleContr[A_contrSet])
          if(length(B_dropoutSet)>0) likEvid = likEvid*prod(dropAlleleContr[B_dropoutSet])
          if(length(C_dropinSet)>0) {
            likEvid = likEvid*prod(pC*freq[C_dropinSet])
          } else {
            likEvid = likEvid*(1-pC) #no dropin
          }
        }
        #FINALLY INCLUDE PERMUTATION FACTOR AND UPDATE LIKELIHOOD
        genoComb =  exp(m_genoCombConstant - sum(lgamma( table(genoJointIdx)+1 )))
        
        liksum <<- liksum + likEvid*genoProb*genoComb #genoComb 
        
      } else {
        recfun(contrIdx+1, genoJointIdx, genoProb, contrAlleles, nTyped)
      }
      
    }  
    return()
  } 

  #tmp = lgamma(m_nGenos + nUnknowns - 1 + 1) - lgamma(m_nGenos - 1 + 1) - lgamma(nUnknowns + 1)
  #nGenosJointRestricted = as.integer(exp(tmp))
  liksum <<- 0 #BIG SUM
  if(nUnknowns>0) {
    recfun(1, rep(0,nUnknowns), 1, m_contrAlleles, m_nTyped)
  } else {
    hasContr = m_contrAlleles > 0 #bool of which has contributin
    hasNoContr = !hasContr
    dropAlleleContr = pD^m_contrAlleles#[hasContr] #calculate log  pD^n_a (used two places)
    
    likEvid = 1
    for(r in seq_len(nReplicates)) {
      hasEvid = hasEvidMat[,r]
      A_contrSet = which(hasEvid & hasContr) #Contribution set
      B_dropoutSet = which(!hasEvid & hasContr) #Dropout set
      C_dropinSet = which(hasEvid & hasNoContr) #Dropin set
      if(length(A_contrSet)>0) likEvid = likEvid*prod(1-dropAlleleContr[A_contrSet])
      if(length(B_dropoutSet)>0) likEvid = likEvid*prod(dropAlleleContr[B_dropoutSet])
      if(length(C_dropinSet)>0) {
        likEvid = likEvid*prod(pC*freq[C_dropinSet])
      } else {
        likEvid = likEvid*(1-pC) #no dropin
      }
    }
    liksum=likEvid
  }
  
  lik1 = as.numeric(liksum)
  return(lik1)
  #nContr = length(conds)/2 + nUnknowns
  #lik2 = forensim::likEvid(evids,conds,knowns,nUnknowns,fst, rep(pD,nContr),rep(pD,nContr)^2,pC,freq)
  
 }
  
  
  