#This function contains a bag of helpfunctions used for some functions in the package.


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
      if(!isPhi) param2[,idx] = phi[,NOC+cc] #add if found
      if(isPhi) param2[,idx] = .invlogit(phi[,NOC+cc]) #add if found
      cc = cc + 1
    }     
  }
  param1 = phi[,c(NOC,NOC+1)] #obtain PHexp/PHvar
  if(isPhi) param1 = exp(param1) #transform back
  
  #Output depends on dimension
  if(dim==1) {
    return( c( mixprop,param1,param2) )
  } else {
    return( cbind( mixprop,param1,param2[,modTypes,drop=FALSE]) ) #last params are dropped (not used in MCMC output)
  }
}

#Helpfunction to draw "good" randoms start points (good guess)
.paramrandomizer = function(th0,NOC,modTypes,delta,ncond=0,hasKinship=0,verbose=FALSE) { 
  .logit = function(x) log(x/(1-x))
  mxrnd = rgamma(NOC,1) #Draw simplex (flat)
  mxrnd = mxrnd/sum(mxrnd)
  nU = NOC-ncond
  if( nU>1) { #sort if more than 1 unknown
    indSort = (ncond+1):NOC #indicate indices to sort Mix-prop for the unknowns 
    if(as.integer(hasKinship)>0) { #we need to shorten which indices that are sorted
      rmRange = length(indSort) #remove last contributor (this is an unknown related)
      if(nU > 2 && as.integer(hasKinship)==2) rmRange = rmRange - 1:0 #also remove second last (this is also unrelated if hasKinship=2)  
      indSort = indSort[-rmRange] #remove index of last unknown (last element)
    }
    if(length(indSort)>1) mxrnd[indSort] = sort(mxrnd[indSort],decreasing = TRUE) #sort Mx in decreasing order if at least 2
  }

  #convert Mx values to real domain (nu:
  nurnd = numeric()
  if(NOC>1) {
    cs = 0 #c( 0,cumsum(mxrnd)) #cumulative sum of mixture proportins
    for(cc in 1:(NOC-1)) { #traverse contributors (Restricted)
      nurnd = c(nurnd, .logit( mxrnd[cc]/(1-cs))) 
      cs = cs + mxrnd[cc] #update sum
    }
  }
  
  #CONSIDER PH prop variables
  th1 = th0[1:2]  #PH prop variables
  sdPH = delta*0.15*th1 #obtain considered SD of PH props 
  PHrnd = abs( rnorm(2,th1,sd=sdPH))
  logPHrnd = log(PHrnd)  #Obtain random start for mu/sigma, Note using the delta here (should be small)
  randParam = c(mxrnd,PHrnd) #add random for PHprop
  
  #CONSIDER other variables (degrad/BW/FW)
  otherrnd = numeric() #random for beta,xi etc
  if(modTypes[1]) {
    maxval = 5 #maximum transformed degrad slope param
    degval = .logit( th0[3] ) #get transformed degrad slope param value
    if( is.infinite(degval) || degval>maxval ) degval = maxval #insert fixed large value (maxval=5)
    otherrnd = c(otherrnd, degval ) #extract for degradation
  }
  if(modTypes[2]) otherrnd = c(otherrnd, .logit(0.05) ) #assume stutter prop expectation of 0.05
  if(modTypes[3]) otherrnd = c(otherrnd, .logit(0.01) ) #assume stutter prop expectation of 0.025
  if(length(otherrnd)>0) otherrnd = rnorm(length(otherrnd),otherrnd,sd=0.5) #Note small sd (because shifted with expected trend)
  randParam = c(randParam, 1/(1+exp(-otherrnd))) #transform back and add
  if(verbose)   cat(paste0("\nRandom startparam=",paste0(prettyNum(randParam),collapse="/")))
  return(c(nurnd,logPHrnd,otherrnd))    
}

#Delta-method to Sigma matrix
#phi=maxPhi;theta=thetahat
.calcJacobian <- function(phi,theta,NOC) { #Jabobian matrix (in value phi)
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
  
  if(  length(otherInd)>0   ) { #if remaining params: 
    tmp = exp(-phi[otherInd])
    J[cbind(otherInd,otherInd)] = tmp*(1+tmp)^(-2) #insert values (vectorized)
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


