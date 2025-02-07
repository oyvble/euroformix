#A bag of helpfunctions used for some functions in the package

.getThetaUnknowns = function(thetaFull,NOC,modTypes,inclLastMx=FALSE) {
  thetaStartUnknowns = thetaFull[1:(NOC+2)] #obtain Mx and PH vars
  if(!inclLastMx) thetaStartUnknowns = thetaStartUnknowns[-NOC] #remove last Mx
  for(modTypeIdx in seq_along(modTypes)) {
    if(modTypes[modTypeIdx]) thetaStartUnknowns = c(thetaStartUnknowns,thetaFull[NOC+2+modTypeIdx])
  }
  return(thetaStartUnknowns)  
}

#nKnown=0;gridSize=5;epsilon=0.01;last=TRUE
#Helpfunction to get an outcome of mixture proportion set (depends on gridsize)
#nUnknown=1;nKnown=3;gridSize=5;epsilon=0.01;last=TRUE
.getMxOutcome = function(nUnknown,nKnown=0,gridSize=5,epsilon=0.01,last=TRUE) {
  #  gridSize=3
  nTotal =nUnknown+nKnown #get total number of contributors 
  if(nTotal==0) return(NULL)
  
  prepareMxVec = function(x) {
    x[1] = epsilon
    x[length(x)] = x[length(x)] - epsilon
    return(x)
  }
  
  #Obtain list of Mx vectors for each contributor 
  mxList = list() 
  if(nTotal==1) return(matrix(1,nrow=1,ncol=1)) #return empty if one contributor (mx=1 always)
  if(nUnknown>=2) { #if at least two unknowns
    for(k in 2:nUnknown) {
      mxVec = prepareMxVec(seq(0,1/k,l=gridSize)) #construct mx.
      mxList[[k-1]] = mxVec
    }
  }
  if(nKnown>0) {
    mxListKnown = list()
    mxVec = prepareMxVec(seq(0,1,l=2*gridSize)) #construct mx.
    for(k in seq_len(nKnown)) {
      mxListKnown[[k]] = mxVec
    }
    mxList = c(mxListKnown,mxList)
    if(nUnknown==0) mxList = mxList[-1] #remove one of the contributors
  }
  
  #Obtain the matrix outcome
  mxOut = expand.grid(mxList) #obtain outcome
  if(nUnknown>=3) { #if at least 3 unknowns (can restrict)
    urange = nKnown + 2:(nUnknown-1) 
    for(k in urange) {
      keep = mxOut[,k-1]>=mxOut[,k]
      mxOut = mxOut[keep,,drop=FALSE]
    }
  }
  lastMx = 1-rowSums(mxOut) #last is the unknown that is constrained
  
  #Block to avoid invalid restricted Mx
  isValid = lastMx>=0 & lastMx <= 1
  mxOut = mxOut[isValid,,drop=FALSE]
  lastMx = 1-rowSums(mxOut) #last is the unknown that is constrained
  
  if(nUnknown>=2) { #if at least 2 unknowns we want the last unknown to be greatest
    keep = lastMx >= mxOut[,nKnown+1]
    mxOut = mxOut[keep,,drop=FALSE]
  } 

  #Need to make sure that last Mx is not too small
  lastMx = 1-rowSums(mxOut) #obtain Mx of last proportion
  isSmaller = lastMx<epsilon #whether it was smaller
  if(any(isSmaller)) {
    mxOut[isSmaller,] = mxOut[isSmaller,]*(1-epsilon) #updated
  }
  
  if(last) { #if adding last Mx as well
    mxOut = cbind(mxOut, 1-rowSums(mxOut) ) #add last contributor
  }
  rownames(mxOut) = NULL
  colnames(mxOut) = NULL
  return(mxOut)
}

#Performs presearch of the MixProp set (used before calling logLik in MLE, INT and MCMC)
#BWstutt0=0.06; FWstutt0=0.02; epsilon_mx=0.01; epsilon_xi=0.001 
.preSearch = function(obj,c,thetaPresearch,resttol,priorBWS, priorFWS, BWstutt0=0.06, FWstutt0=0.02, epsilon_mx=0.01, epsilon_xi=0.001 ) {
  modTypes = c$modTypes
  
  precalcLogLik = function(param) {
    loglik = obj$calcGenoWeightsMax(as.numeric(param))
    if(modTypes[2] && !is.null(priorBWS))  loglik = loglik + log( priorBWS(param[c$nC+4]) ) #weight with stutter prrior
    if(modTypes[3] && !is.null(priorFWS))  loglik = loglik + log( priorFWS(param[c$nC+5]) ) #weight with stutter prrior 
    return(loglik)
  }    
  
  #CREATE PARAMETER OUTCOME TO CHECK
  nKnown = c$nC - c$nU #number of conditionals
  nUnknown = c$nU #number of Mx params to be restricted
  nKinship = as.integer(c$hasKinship) #number of kinship specifications
  if(nKinship>0) { #must re-adjust number of unknown/known since kinships must vary freely
    nUnknown = nUnknown - nKinship #reduce by one (or two)
    nKnown = nKnown + nKinship
  } 
  
  #Obtain mx outcome
#  .getMxOutcome(nUnknown=1,nKnown=3)
  mxOutcomeSearch = .getMxOutcome(nUnknown,nKnown,5,epsilon_mx)
  numChangeStuttType = 1 + sum(as.integer(modTypes[2:3]))
  presearchFitList = list() #storage for prefit params
  for(changeStuttTypeIdx in 1:numChangeStuttType) {
#    changeStuttTypeIdx=1
    presearchFit = NULL #storage for prefit params
    if(changeStuttTypeIdx==2) {
      #        if(verbose) print(paste0("...pre-searching with other stutter params..."))
      loglik0vec = presearchFitList[[1]][,ncol(presearchFitList[[1]])]
      lik0vec = exp(loglik0vec-max(loglik0vec)) #avoid overflow by subtract max
      lik0vec = lik0vec/sum(lik0vec)
      ord = order(lik0vec,decreasing = TRUE)
      lik0vecOrdered = lik0vec[ord]
#      barplot(lik0vecOrdered)
      ordKeepIdx = cumsum(lik0vecOrdered)<0.9999#set threshold
      mxOutcomeSearch = mxOutcomeSearch[ord[ordKeepIdx],,drop=FALSE]
    }
    for(mxIter in seq_len(nrow(mxOutcomeSearch))) {
      #if(verbose) print(paste0("",round(mxIter/mxSearchOutcome*100),"%.."))
      #        mxIter =    1    
      mxFull = unlist(mxOutcomeSearch[mxIter,]) #already full vector
      #mxRestricted = 1-sum(mx) #this is the restricted mix prop, here defined as the greatest 
      #mxFull = c(mx,mxRestricted)
      
      #Need to special handle situation with kinship (related are the last unknown(s))
      if(nKinship>0) {
        kinshipRange = 1:nKinship #obtain mx range of kinship individuals
        mxFull = c(mxFull[-kinshipRange],mxFull[kinshipRange]) #put these last instead of first (IMPORTANT!)
      }
      param0 = c(mxFull,thetaPresearch) #include restricted mixprop and prefitted params
      if(changeStuttTypeIdx>1) { #reduce stutter prop values if having stutter models
        if(changeStuttTypeIdx<numChangeStuttType) {
          FWstutt0 = epsilon_xi #set as low value
        } else {
          FWstutt0 <- BWstutt0 <- epsilon_xi #set as low value
        }
      }
      param0 = c(param0, ifelse(modTypes[2], BWstutt0,0)) #always append 
      param0 = c(param0, ifelse(modTypes[3], FWstutt0,0)) #always append 
      #param0[1:2] = 0.39; param0[3:4] = 0.11
      #param0[1:4]=param0[1:4]/sum(param0[1:4])
      logLik0 = precalcLogLik(param0) #NOTE: Unstability if param not given as full vector (zeros must be added!)
      presearchFit = rbind(presearchFit, c(param0, logLik=logLik0)) #search object
    }
    #plot(as.data.frame(presearchFit[,c(2:nC)]))
    presearchFitList[[changeStuttTypeIdx]] = presearchFit
  } #end for each stutter model type
  
  #Step 4: Perform restriction and give update:
  resttolReps = resttol^c$nReps #make more restriction when replicates (squared)
  restprop = obj$restrictGenos(as.numeric(resttolReps))  
  return(list(presearchFitList=presearchFitList,restprop=restprop))
}  #end function



#Helpfunctions for the transformation between the theta and phi domain for deg/stutt params
#.theta2phi = function(x) log(1+exp(x)) #softplus
#.phi2theta = function(x) log(exp(x)-1) #softplus inverted
.theta2phi = function(x) log(x)
.phi2theta = function(x) exp(x)

.getMxValid = function(mxVec,nC,nU,hasKinship) {
  nKinship = as.integer(hasKinship)
  nUnknown = nU - nKinship #count only unrelated unknowns
  nCond = nC - nU #number of conditionals (not kinships)
  #Note that Kinship unknowns are expected last (can vary free)
  mxVecValid = mxVec
  if(nUnknown>=2) { #if at least 2 unknowns unrelated (can restrict)
    urange = nCond + 1:nUnknown #range of these
    mxVecU = sort(mxVec[urange],decreasing = TRUE) #make sort of unknown part
    mxVecValid[urange] = c(mxVecU[-1],mxVecU[1]) #ensure largest last
  }   
  return(mxVecValid) #return potential modified copy
}

#helpfunction to check if the order of the Mx vector is OK:
.checkMxValid = function(mxVec,nC,nU,hasKinship) {
  nKinship = as.integer(hasKinship)
  nUnknown = nU - nKinship #count only unrelated unknowns
  nCond = nC - nU #number of conditionals (not kinships)
  #Note that Kinship unknowns are expected last (can vary free)
  
  mxValid = TRUE
  if(nUnknown>=2) { #if at least 2 unknowns unrelated (can restrict)
    urange = nCond + 1:nUnknown #range of these
    mxU1 = mxVec[urange[1]] #get first U
    mxUL = mxVec[urange[length(urange)]] #get last U
    if( mxUL < mxU1) { #last Mx of unknown cannot be smaller than first
      mxValid = FALSE
    } else if(nUnknown>=3) { #otherwise another situation with at least 3 unknowns
      urangeNotFirstLast = urange[-c(1,length(urange))] #remove first and last element
      for(k in urangeNotFirstLast) {
        if( mxVec[k] > mxVec[k-1] ) { #check if current is greater than prevous
          mxValid = FALSE
          break
        }
      }
    }
  }
  return(mxValid)
}


#Helpfunction in calcMLE to get covariances and names
.getFittedParams = function(maxPhi,maxSigma,NOC,modTypes) {
  np = length(maxPhi)
  theta2full = .convBack(maxPhi,NOC,modTypes) #include restricted 
  theta2 = theta2full[1:(NOC+2)]
  for(modTypeIdx in seq_along(modTypes)) {
    if(modTypes[modTypeIdx]) theta2 = c(theta2,theta2full[NOC+2+modTypeIdx])
  }
  theta = theta2[-NOC]
  
  J <- .calcJacobian(maxPhi,theta,NOC) #obtain jacobian matrix:
  Sigma <- (t(J)%*%maxSigma%*%J) #this is correct covariance of thetahat. Observed hessian is used
  
  #OBTAIN COVARIANCE MATRIX OF PARAMS:
  #get extended Sigma (all parameters)
  Sigma2 <- matrix(NA,nrow=np+1,ncol=np+1) #extended covariance matrix also including mx[nC]
  Sigma2[NOC:np+1,NOC:np+1] <- Sigma[NOC:np,NOC:np] 
  Sigma2[NOC:np+1,NOC:np+1] <- Sigma[NOC:np,NOC:np] 
  if(NOC>1) {
    Sigma2[NOC:np+1,1:(NOC-1)] <- Sigma[NOC:np,1:(NOC-1)] 
    Sigma2[1:(NOC-1),1:(NOC-1)] <- Sigma[1:(NOC-1),1:(NOC-1)] 
    Sigma2[1:(NOC-1),NOC:np+1] <- Sigma[1:(NOC-1),NOC:np] 
    Sigma2[NOC,NOC] <- sum(Sigma[1:(NOC-1),1:(NOC-1)])
    for(k in (1:(np+1))[-NOC]) {
      Sigma2[NOC,k] <- Sigma2[k,NOC] <- -sum(Sigma[1:(NOC-1),k-sum(k>NOC)]) 
    }
  } else {
    Sigma2[1,1:(np+1)] <- Sigma2[1:(np+1),1] <- 0 #no uncertainty
  }
  
  #Standard error for theta:
  suppressWarnings({
    thetaSE <- sqrt(diag(Sigma2))
    thetaSE[is.nan(thetaSE)] = 0 #avoid NAN
  })
  
  mxName <- "Mix-prop. C" #"mx"  
  thetanames0 <- c("P.H.expectation","P.H.variability")
  phinames <- paste0("log(",thetanames0,")")
  if(NOC>1) {
    phinames  <- c(paste0("nu",1:(NOC-1)),phinames)
    thetanames <- c(paste0(mxName,1:(NOC-1)),thetanames0)
  } else {
    thetanames <- thetanames0
  }
  thetanames2 <- c(paste0(mxName,1:NOC),thetanames0)
  
  othernames <- character() #also include stutter models
  if(modTypes[1]) othernames <- c(othernames,"Degrad. slope")
  if(modTypes[2]) othernames = c(othernames,"BWstutt-prop.")  
  if(modTypes[3]) othernames = c(othernames,"FWstutt-prop.")
  if(length(othernames)>0) phinames <- c(phinames,paste0("log(",othernames,")") )
  thetanames <- c(thetanames, othernames )
  thetanames2 <- c(thetanames2, othernames )
  
  colnames(maxSigma) <- rownames(maxSigma) <- phinames 
  colnames(Sigma) <- rownames(Sigma) <- thetanames
  colnames(Sigma2) <- rownames(Sigma2) <- thetanames2
  names(maxPhi) <- phinames
  names(theta) <- thetanames
  names(theta2) <- thetanames2
  names(thetaSE) <- thetanames2
  return(list(phihat=maxPhi,thetahat=theta,thetahat2=theta2,phiSigma=maxSigma,thetaSigma=Sigma,thetaSigma2=Sigma2,thetaSE=thetaSE))
}


#helpfunction to transfer from theta to phi
.getPhi = function(theta,NOC,modTypes) {
  .logit = function(x) log(x/(1-x))
  mxvec = theta[1:(NOC-1)]
  nuvec = numeric()
  if(NOC>1) {
    cs = 0 #c( 0,cumsum(mxrnd)) #cumulative sum of mixture proportins
    for(cc in 1:(NOC-1)) { #traverse contributors (Restricted)
      nuvec = c(nuvec, .logit(mxvec[cc]/(1-cs))) 
      cs = cs + mxvec[cc] #update sum
    }
  }
  PHvars = log(c(theta[NOC + 0:1]))
  othervars = numeric() #random for beta,xi etc
  if(modTypes[1]) othervars = c(othervars, .theta2phi(theta[NOC + 2]) ) #extract for degradation
  if(modTypes[2]) othervars = c(othervars, .theta2phi(theta[NOC + 2 + sum(modTypes[1])]) ) #assume stutter prop expectation of 0.05
  if(modTypes[3]) othervars = c(othervars, .theta2phi(theta[NOC + 2 + sum(modTypes[1:2])]) ) #assume stutter prop expectation of 0.025
  return(c(nuvec,PHvars,othervars))
}


#HELPFUNCTIONS FOR TRANSFER VARIABLES Between Real<->[0,1] domains
.convBack = function(phi,NOC,modTypes, isPhi=TRUE) {
  .invlogit = function(x) 1/(1+exp(-x))
  dim = nrow(phi)
  if(is.null(dim)) {
    dim = 1
    phi = rbind(phi)
  }  
  mixprop = matrix(1,nrow=dim,ncol=NOC)
  if(NOC>1) {
    cs = 0 #cumulative summing of mixture proportions
    for(i in 1:(NOC-1)) { #for each mix props
      if(!isPhi) mixprop[,i] = phi[,i] #no transform (direct include)
      
      if(isPhi) {
        mixprop[,i] = .invlogit(phi[,i]) #transform from R to M
        if(i>0) {
          mixprop[,i] = mixprop[,i]*(1-cs); #transform further (see formula for explanation)
        }
      }
      cs = cs + mixprop[,i] #add mixture proportion
    }
    mixprop[,NOC] = 1-cs
  } #end if more than 1 contr
  
  #Prepare param2 first
  param2 = cbind(1, matrix(0,nrow=dim,ncol=2)) #DEG,BWS,FWS
  cc = 2 #counter index for additional params
  for(idx in 1:3) { #looping all types
    if(modTypes[idx]) {
      if(isPhi) {
        param2[,idx] = .phi2theta(phi[,NOC+cc]) #add if found
      } else {
        param2[,idx] = phi[,NOC+cc] #add if found
      }
      cc = cc + 1
    }     
  }
  param1 = phi[,c(NOC,NOC+1)] #obtain PHexp/PHvar
  if(isPhi) param1 = exp(param1) #transform back
  
  #Output depends on dimension
  if(dim==1) {
    return( c( mixprop,param1,param2) )
  } else {
    #For vectorized version we will only return unknowns (MCMC)
    return( cbind( mixprop,param1,param2[,modTypes,drop=FALSE]) ) #last params are dropped (not used in MCMC output)
  }
}

#Delta-method to Sigma matrix
#phi=maxPhi;theta=thetahat
.calcJacobian <- function(phi,theta,NOC,hasDEG) { #Jabobian matrix (in value phi)
  otherInd = numeric()
  if(length(phi) > (NOC+1)) otherInd = (NOC+2):length(phi) #obtain index of more indices
  
  J <- matrix(0,length(phi),length(phi))
  if(NOC>1) {
    DmDm <- matrix(0,NOC-1,NOC-1)
    irange <- 1:(NOC-1) #range of mixture proportions
    xtmp <- 1/(1+exp(-phi[irange ])) #aux variable
    dxtmp <-  exp(-phi[irange])*xtmp^2 #derivative of aux variable
    for(i in irange) {
      for(j in 1:i) {
        if(j==i) {
          DmDm[i,i] <- dxtmp[i]
          if(i>1) DmDm[i,i] <- DmDm[i,i]*(1-sum(theta[1:(j-1)])) #note using mle(theta) here!
        } else { #cross-derivatives
          DmDm[i,j] <-  DmDm[j,i] <- -xtmp[i]*sum(DmDm[1:(i-1),j])
        }
      } #end for each col j (xj)
    } #end for each row i (fi)
    J[1:(NOC-1),1:(NOC-1)] <- DmDm
  }
  for(i in NOC:(NOC+1)) J[i,i] <- exp(phi[i])
  
  if(length(otherInd)>0) { #if remaining params: 
    J[cbind(otherInd,otherInd)] = exp(phi[otherInd])
  }
  return(J)
} #end jacobian

.secondToTimeformat <- function(t){ #converts seconds to time format
  paste(formatC(t %/% (60*60) %% 24, width = 2, format = "d", flag = "0"), #hours
        formatC(t %/% 60 %% 60, width = 2, format = "d", flag = "0"), #mins
        formatC(t %% 60, width = 2, format = "d", flag = "0"),sep = ":") #seconds
}


#Helpfunction to print tables (from R v3.5.0) has problem showing tables in the GUI.
.printTable = function(x) {
  print(cbind(rownames(x),x),row.names=FALSE)
}


#helpfunction to get small number from log-value
.getSmallNumber = function(logval,sig0=2,scientific="e") {
  if(is.nan(logval)) {
    return(NaN) #return NAN if not able to calculate
  } else if(is.infinite(logval)) {
    return(0)
  } else if(round(logval,sig0)==0) {
    return(1)
  }
  log10base = logval/log(10) #convert to 10 base
  power = floor(log10base) #get power number
  remainder = log10base - power
  return( paste0( round(10^remainder,sig0),scientific,power)) #representation of very small numbers (avoid underflow)
}

#Helpfunction to ensure correct structure of refData: refData[[rr]][[loc]]
#used in plotEPG/plotTopEPG funcs
.getRefData = function(refData=NULL,locs) { #
  refData2 = NULL
  if(!is.null(refData)) {
    refData2 = list()
    if(any(names(refData)%in%locs)) { #assuming this structure refData[[loc]][[rr]] 
      refNames = unique(unlist(lapply(refData,names))) #get reference names
      for(ref in refNames) refData2[[ref]] = lapply(refData,function(x) x[[ref]])  #converting structure
      #Note that refData structure may have been changed if missing markers (irregularity)
    } else {
      refData2 = refData #format was OK
    }
  }
  return(refData2)
}

#Helpfunction to advice number of MCMC iter (from mcmcse)
#https://github.com/cran/mcmcse/blob/master/R/minESS.R
.getMinESS = function(p, alpha = .05, eps = .05) {
  crit <- qchisq(1-alpha, p)
  foo <- 2/p
  logminESS <- foo*log(2) + log(pi) - foo*log(p) - foo*lgamma(p/2) - 2*log(eps) + log(crit)
  return(round(exp(logminESS)))
}

.getPrim = function() as.integer(c(2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,103,107,109,113, 127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223,227,229,233,239,241,251,257,263, 269,271,277,281,283,293,307,311,313,317,331,337,347,349,353,359,367,373,379,383,389,397,401,409,419,421, 431,433,439,443,449,457,461,463,467,479,487,491,499,503,509,521,523,541,547,557,563,569,571,577,587,593, 599,601,607,613,617,619,631,641,643,647,653,659,661,673,677,683,691,701,709,719,727,733,739,743,751,757, 761,769,773,787,797,809,811,821,823,827,829,839,853,857,859,863,877,881,883,887,907,911,919,929,937,941, 947,953,967,971,977,983,991,997,1009,1013,1019,1021,1031,1033,1039,1049,1051,1061,1063,1069,1087,1091,1093, 1097,1103,1109,1117,1123,1129,1151,1153,1163,1171,1181,1187,1193,1201,1213,1217,1223,1229,1231,1237,1249, 1259,1277,1279,1283,1289,1291,1297,1301,1303,1307,1319,1321,1327,1361,1367,1373,1381,1399,1409,1423,1427, 1429,1433,1439,1447,1451,1453,1459,1471,1481,1483,1487,1489,1493,1499,1511,1523,1531,1543,1549) )

#helpfunction to get transparant contribution colors
.getContrCols <- function(deg = .3) { 
  colNames <- c("blue","green","red","yellow","pink","brown","orange","lightgoldenrod") #contributor cols
  colNamesT = setNames(colNames,colNames)
  for(color in colNames) {
    color2 = color
    switch(color,
           "green" = {color2 = "green4"},
           "red" = {color2 = "coral3"},
           "orange" = {color2 = "darkorange"},
           "brown" = {color2 = "saddlebrown"},
           "pink" = {color2 = "magenta"},
           "lightgoldenrod" = {color2 = "lightgoldenrod3"},
           "yellow" = {color2 = "gold2"})
    rgb.val <- col2rgb(color2)
    colNamesT[color] =  rgb(rgb.val[1], rgb.val[2], rgb.val[3],maxColorValue = 255, alpha = (1-deg)*255)
  }
  return(colNamesT)
}

#Helpfunction to obtain fragment length for given allele
.getFragLength = function(alleles,kittab, isSTRING=FALSE) {
  #kittab is locus specific table
  #alleles <<- alleles
  #kittab <<- kittab
  #isSTRING <<- isSTRING
  
  fragLength = kittab$Size[match(alleles,kittab$Allele)] #corresponding fragment length of alleles
  names(fragLength) = alleles
  idxNAfragLength = which(is.na(fragLength)) #obtain indices of alleles that did not find fragmentlenght
  
  if(length(idxNAfragLength)>0) {
    if(isSTRING) { #if string
      for(i in idxNAfragLength) fragLength[i] =  mean(kittab$Size) #take simple average
    } else {
      #fit regression line:
      isWholeNumber = !grepl("\\.",kittab$Allele)
      if(sum(isWholeNumber)<2) { #if less than two
        for(i in idxNAfragLength) fragLength[i] =  kittab$Size[which.min(abs(as.numeric(kittab$Allele) - as.numeric(alleles[i])))]
      } else {
        #plot(alleles_num[!isNAfragLength],fragLength[!isNAfragLength])
        CEname = as.numeric(kittab$Allele[isWholeNumber][1:2])
        CEbp = as.numeric(kittab$Size[isWholeNumber][1:2])
        motiflen = round(diff(CEname)*diff(CEbp)) #obtain motif length
        offset = CEbp[1] - (CEname[1]*motiflen)
        
        predAlleleList = strsplit(as.character(alleles[idxNAfragLength]),"\\.")
        for(i in seq_along(predAlleleList)) {
          bp = as.numeric(predAlleleList[[i]][1])*motiflen #obtain fragment length
          if(length(predAlleleList[[i]])>1) bp = bp + as.numeric(predAlleleList[[i]][2]) #add number
          fragLength[idxNAfragLength[i]] = bp + offset #insert final
        }
      }
    }
  }
  #plot(as.numeric(names(fragLength)),fragLength)
  return(fragLength)
}


#helpfunction to obtain data for showing in plotTop genotypes
.getDataToPlotProfile = function(mlefit, DCobj, kit=NULL, withStutterModel=TRUE, grpsymbol=NULL,Qallele = "99") {
  c = mlefit$prepareC #obtain C-preparation object
  
  #extract info from DC (deconvolution) object
  if(is.null(DCobj)) DCobj <- deconvolve(mlefit,maxlist=1) #get top candidate profiles
  if(is.null(DCobj)) return(NULL) #could not be done
  topGlist <- DCobj$toprankGi #get top genotype info
  
  #obtain estimates:
  model = mlefit$model
  AT = model$threshT #extract analytical threholds from object
  thhat <- mlefit$fit$thetahat2 #get estimates
  nC <- c$nC
  mx <- thhat[1:nC]   
  mu <- thhat[nC+1]   
  sigma <- thhat[nC+2]   
  beta <- 1
  xiBW <- xiFW <- 0
  if(c$modTypes[1]) beta <- thhat[nC+3] #DEG
  if(c$modTypes[2]) xiBW <- thhat[nC+3 + sum(c$modTypes[1]) ] #BW stutter
  if(c$modTypes[3]) xiFW <- thhat[nC+3 + sum(c$modTypes[1:2]) ] #FW stutter
  
  nReps = c$nReps #number of replicates
  locsAll = c$markerNames #get locus names
  sampleNames = c$repNames #obtain evidence names
  refNames = c$refNamesCond #obtain reference name of conditionals
  nrefs = length(refNames)
  #nrefs==c$NOK
  
  #Create dataset (per dye info with bp)
  df = NULL #store data: (sample,marker,allele,height,bp)
  for(locidx in seq_along(locsAll)) {
    #  locidx = 1
    loc = locsAll[locidx]
    
    #Obtain allele/peak/fragmentLen details for marker:
    genoNames = c$genoList[[loc]]
    
    nAlleles = c$nAlleles[locidx] #number of alleles (not replicates)
    nReps = c$nRepMarkers[locidx]
    
    offset_alleles = c$startIndMarker_nAlleles[locidx]
    offset_allelesReps = c$startIndMarker_nAllelesReps[locidx]
    
    alleles = c$alleleNames[offset_alleles + seq_len(nAlleles)] #obtain allele names
    peaksAlleles = c$peaks[offset_allelesReps + seq_len(nAlleles*nReps)]
    peaksAlleles = matrix(peaksAlleles,nrow=nReps,ncol=nAlleles)
    bpAlleles = c$basepairs[offset_alleles + seq_len(nAlleles)]
    
    #Obtain genotypes to visualize
    knownRefGeno = matrix("",ncol=2,nrow=nrefs)
    if(nrefs>0) {
      ginds = c$knownGind[nrefs*(locidx-1) + seq_len(nrefs)] + 1 #note adjust index
      ins = which(ginds>0)
      if(length(ins)>0) knownRefGeno[ins,] = genoNames[ginds[ins],,drop=FALSE] #insert if non-empty
    }
    
    #ref text under each allele
    reftxt = rep("",length(alleles))
    for(rr in seq_len(nrefs)) { #for each ref
      indadd = which(alleles%in%knownRefGeno[rr,]) #index of alleles to add to text
      hasprevval = indadd[nchar(reftxt[indadd])>0] #indice to add backslash (sharing alleles)
      reftxt[ hasprevval ] = paste0(reftxt[ hasprevval ],"/")      
      reftxt[indadd] = paste0( reftxt[indadd], rr)
    }
    
    #GET cumulative model information, EXPECTation (and std) of PH
    Gcontr = topGlist[[loc]][1,] #get top genotypes
    EXPmat <- matrix(0,nrow=length(alleles),ncol=nC)
    SHAPEvec <- rep(0,length(alleles)) 
    for (aa in seq_along(alleles)) { # Loop over all alleles in locus
      allele = alleles[aa]
      contr <- sapply(strsplit(Gcontr,"/"),function(x) sum(x%in%allele)) #get contribution
      EXPmat[aa,] <- cumsum(contr*mx*mu*beta^bpAlleles[aa]) #expected peak heights for each contributors
      SHAPEvec[aa] <- sum(contr*mx)*beta^bpAlleles[aa]/sigma^2 #expected peak heights for each contributors
    }
    
    nStutt = c$nStutters[locidx] #obtain number of stutters
    if(withStutterModel && nStutt>0) {
      stuttidx = c$startIndMarker_nStutters[locidx] + seq_len(nStutt) #get idx range
      stuttFrom = c$stuttFromInd[stuttidx] + 1 #adjust index
      stuttTo = c$stuttToInd[stuttidx] + 1  #adjust index
      #stuttParamIdx = c$stuttParamInd[stuttidx]+1 #obtain parameter indx
      stuttProp = c(xiBW,xiFW)[c$stuttParamInd[stuttidx]+1] #obtain stutter proportions
      
      #scale expectation and shape with stutter:
      EXPmat2 = EXPmat
      SHAPEvec2 = SHAPEvec
      for(ss in seq_len(nStutt) ) {
        from = stuttFrom[ss]
        to = stuttTo[ss]
        EXPmat2[from,] = EXPmat2[from,] - stuttProp[ss]*EXPmat[from,] #SUBTRACTED stutters
        SHAPEvec2[from] = SHAPEvec2[from] - stuttProp[ss]*SHAPEvec[from]  
        if( to <= length(alleles)) {
          EXPmat2[to,] = EXPmat2[to,] + stuttProp[ss]*EXPmat[from,] #OBTAINED stutters
          SHAPEvec2[to] = SHAPEvec2[to] + stuttProp[ss]*SHAPEvec[from]  
        }   
      }
      EXPmat = EXPmat2 #override with stutter products
      SHAPEvec = SHAPEvec2
    }
    PIlow = qgamma(0.025, shape=SHAPEvec, scale=mu*sigma^2)
    PIup = qgamma(0.975, shape=SHAPEvec, scale=mu*sigma^2)
    EXP <- apply(EXPmat,1,function(x) paste0(x,collapse = "/"))
    
    #SPECIAL HANDLING FOR MPS (use group symbol to extract)
    allelesCE = alleles
    if(!is.null(grpsymbol)) {
      splitlist = strsplit(alleles,grpsymbol)  #split allele wrt split symbol
      allelesCE = sapply(splitlist,function(x) x[1]) #extract RU allele
    }
    
    #Store values in table (not Q-allele if not relevant)
    bpAlleles2 = (100*bpAlleles) + 125 #adjust back
    Qidx = which(alleles%in%Qallele) #get index of Q-allele 
    showQ =  reftxt[Qidx]!="" ||  SHAPEvec[Qidx]>0 #whether Q-allele should be included
    if(!is.na(showQ) && !showQ) {
      alleles = alleles[-Qidx]
      allelesCE = allelesCE[-Qidx]
      peaksAlleles = peaksAlleles[,-Qidx,drop=FALSE]
      bpAlleles2 = bpAlleles2[-Qidx] 
      reftxt = reftxt[-Qidx] 
      EXP = EXP[-Qidx] 
      PIlow = PIlow[-Qidx]
      PIup = PIup[-Qidx]
    }
    #Obtaining genotypeProbabilties (same for whole marker)
    genoProb = min(as.numeric(topGlist[[loc]][2,])) #get smallest prob
    
    for(ss in seq_len(nReps)) {
      df = rbind(df, cbind(sampleNames[ss],loc,alleles,allelesCE,peaksAlleles[ss,],bpAlleles2,reftxt,EXP,PIlow,PIup, genoProb=genoProb) )
    }
  } #end for each loci
  df = data.frame(Sample=df[,1],Marker=df[,2], Allele=df[,3],AlleleCE=as.character(df[,4]),Height=as.numeric(df[,5]),bp=as.numeric(df[,6]),reftxt=df[,7],EXP=df[,8],PIlow=as.numeric(df[,9]),PIup=as.numeric(df[,10]),genoProb=as.numeric(df[,11]),stringsAsFactors=FALSE)
  return(df)
}

