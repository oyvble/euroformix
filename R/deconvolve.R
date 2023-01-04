#' @title deconvolve
#' @author Oyvind Bleka
#' @description deconvolve ranks the set of the most conditional posterior probability of genotypes the STR DNA mixture given a fitted model under a hypothesis.
#' @details The procedure calculates the likelihood for each single locus. Then it combines the most probable genotypes from each loci to produce a ranked list of deconvolved profiles.
#' 
#' @param mlefit Fitted object using contLikMLE function.
#' @param alpha Required sum of the listed posterior probabilities.
#' @param maxlist The ranked deconvolved profile list will not exceed this number (used to avoid endless search).
#' @param signif Number of significant numbers
#' @param checkCalcs Return only calculated log-likelihood per markers (summed). Used for testing.
#' @return ret A list(table1,table2,table3,table4,rankGi,rankG,pG) where rankG is the ranked genotypes with corresponding probabilities in pG. rankgGi is the same, but per marker. table1 is rankG and pG combined (joint results). table2 uses rankGi to find marginal results for top-genotypes. table3 and table4 shows this marginalized on genotypes and alleles per contributor 
#' @export

deconvolve = function(mlefit,alpha=0.95,maxlist=1000,signif=4,checkCalcs=FALSE){ 
  c <- mlefit$prepareC #returned from prepareC
  nM = c$nMarkers #number of markers
  locs = c$markerNames #loci to evaluate
  modelType = mlefit$modelType #obtain model type (DEG,BWS,FWS)
  nC = c$nC #number of contrs
  
  #prepare format parameter (required for other code)
  th = mlefit$fit$thetahat2 #obtain params
  if(any(is.na(th))) return(NULL) #could not deconvolve!
  par = list(mixProp=th[1:nC], PHexp=th[nC+1], PHvar=th[nC+2], DEG=1, stutt=rep(0,2) )
  cc = 3 #counter
  if(modelType[1]) { par$DEG = th[nC+cc]; cc = cc + 1 } #update counter
  if(modelType[2]) { par$stutt[1] = th[nC+cc]; cc = cc + 1 } #update counter
  if(modelType[3]) par$stutt[2] = th[nC+cc]; 
  
  #Step 1) Calculate L(E|g,thetahat) for each marker: THIS IS DONE TROUGH SCRIPT
  #nrow(combGind)==c$nJointCombs[m] #must be the same
  loglikVEC = rep(0,sum(c$nJointGenos)) #init big calculation vector (likelihood for all outcome)
  loglikVEC = .C("loglikGamma_allcomb2", as.numeric(loglikVEC),c$startIndMarker_nJointGenos, c$nJointGenos, c$nC, c$NOK, c$nKnowns,
            as.numeric(par$mixProp),  as.numeric(par$PHexp), as.numeric(par$PHvar), as.numeric(par$DEG), as.numeric(par$stutt), 
            c$AT,c$fst,c$dropinProb,c$dropinWeight, c$nMarkers, c$nRepMarkers, c$nAlleles, c$nPotStutters,
            c$startIndMarker_nAlleles, c$startIndMarker_nAllelesReps, c$peaks, c$freqs, c$nTyped, c$maTyped, c$basepairs, 
            c$nStutters, c$stuttFromInd, c$stuttToInd, c$stuttParamInd , c$startIndMarker_nStutters,
            c$knownGind, c$relGind, c$ibd, as.integer(mlefit$maxThreads) )[[1]] #obtain likelihood of all outcome
  
  #length(loglikVEC)
  #CHECK THAT SAME LIKELIHOODs ARE OBTAINED:
  if(checkCalcs) {
    logLiki = rep(0,nM)
    for(m in 1:nM) logLiki[m] = log(sum(exp( loglikVEC[ c$startIndMarker_nJointGenos[m] + 1:c$nJointGenos[m] ])))
    return(logLiki)
  }

  #Helpfunctions to obtain marginal probabilities
  getMarg <- function(x,y) { #get marginal of genotypes
    agg <- aggregate(y,by=list(x),sum) #get probabilities
    ord <- order(agg[,2],decreasing=TRUE)
    agg2 <- agg[ord,,drop=F]
    colnames(agg2) <- c("Genotype","Probability")
    return(agg2)
  }
  getMarg2 <- function(x,y) { #get marginal of alleles
    tmp <- unlist(strsplit(x,"/"))
    unA <- unique(tmp) #unique alleles
    x2 <- t(matrix(tmp,nrow=2))
    prob <- rep(NA,length(unA))  
    for(aa in unA) prob[which(unA==aa)] <- sum(y[rowSums(x2==aa)>0]) #sum probabilities
    ord <- order(prob,decreasing=TRUE)
    agg <- data.frame(Allele=unA[ord],Probability=prob[ord])
    return(agg)
  }
  
  #helpfunction to obtain a maximum size of a vector (bounded in both length and probability)
  maxI <- function(p) min(min(which(cumsum(p)>=alpha)),maxlist,length(p))   #keeping cumulative sum
  maxI2 <- function(p) min(max(which(p>(1-alpha))),maxlist,length(p)) #keeping absolute sum
  
  
  #Step 2) Obtaining Deconvolved Genotype  results
  contrcn <- paste0("C",1:nC) #column name of contributors
  topRankcn <-  c("TopGenotype","probability","ratioToNextGenotype") #names for each contributor
  
  #OUTPUT VARIABLES
  deconvlisti <- list() #stored full list (truncated on maxlist)
  deconvlistic <- list() #genotype list per contributor
  deconvlistica <- list() #allele list per contributor
  toplist <- list()
  table1 <- table2 <- table3 <- table4 <- numeric()
  for(m in 1:nM) { #extract info in c relevant for each markers:
# m=2;   print(m)
    loc =locs[m] #obtain locus name
    Gset = c$genoList[[loc]]
    Gset1 = paste0(Gset[,1],"/",Gset[,2]) #genotype names
    nGenos1 = length(Gset1)
    NOK0 = c$NOK #obtain number of knowns
    nU = c$nC #number of unknowns
    uind = 1:nC
    kind = integer()
    if(NOK0>0) {
      knownGind = c$knownGind[(m-1)*NOK0+1:NOK0] #obtain genotype index for known contributors
      kind = which(knownGind!= (-1)) #index of knowns
      uind = setdiff(uind,kind) #update
      nU = length(uind) 
    } 
    
    #Obtain outcome of genotypes
    probVec  <- 1 #default probabilty if all are known
    if(nU>0) {
      genoinds <- c$startIndMarker_nJointGenos[m] + 1:c$nJointGenos[m] #get index of genotype outcome
      probVec <- exp(loglikVEC[genoinds]) #obtain calculated values
      probVec <-  probVec/sum(probVec) #normalize
      
      #Obtain genotype combination for unknowns
      unknownGenoJoint = list()
      for(k in seq_len(nU)) unknownGenoJoint[[k]] <- Gset1
      unknownGenoJoint = expand.grid(unknownGenoJoint,stringsAsFactors = FALSE) #get all combinations
    }
    nGenosJoint = length(probVec) #number of joint genotype outcome
    
    #Part 1: Construct joint deconvolved results
    ord = order(probVec,decreasing = TRUE) #Obtain rank of rows
    deconvOBJ <-  matrix(NA,ncol=nC,nrow=nGenosJoint) #translate indices to genotype names
    for(k in 1:nC) { #for each contributors
      if(k%in%kind) genoVec =  rep( Gset1[knownGind[k]+1],nGenosJoint)
      if(k%in%uind) genoVec <- unknownGenoJoint[ord,which(k==uind)] #Obtain unknown genotypes
      deconvOBJ[,k] = genoVec
    }
    probVec = probVec[ord] #sort probabilities
    colnames(deconvOBJ) = contrcn
    
    #Part 2: CALCULATE MARGINAL PRODUCTS (by-products of joint result deconvOBJ)
    #A) Calculate marginal probabilities for all contributors (genotypes and alleles):
    deconvlistica[[loc]] <- deconvlistic[[loc]] <- list()
    tab <- matrix(,nrow=3,ncol=nC)
    rownames(tab) <- topRankcn 
    colnames(tab) <- contrcn
    for(k in 1:nC) {
      deconvlistic[[loc]][[contrcn[k]]] <- tmp <- getMarg(x=deconvOBJ[,k],y=probVec)
      deconvlistica[[loc]][[contrcn[k]]] <- getMarg2(x=deconvOBJ[,k],y=probVec)
      rat <- ifelse(nrow(tmp)>1, tmp[1,2]/tmp[2,2],NA) #get ratio from first to second genotype
      tab[,k] <- c(tmp[1,1],signif(tmp[1,2],signif),signif(rat,signif)) 
    }
    toplist[[loc]] <- tab
    
    #B) Create tables
    maxind <- maxI(probVec) #get max index to use
    if(is.infinite(maxind)) maxind=length(probVec)
    indkeep = seq_len(maxind) #obtain indices to keep
  
    #Obtain limited objects:
    deconvOBJ2 = deconvOBJ[indkeep,,drop=FALSE]
    probVec2 = probVec[indkeep]
    
    combs = cbind(deconvOBJ2,Probability=probVec2) #limited object
    deconvlisti[[loc]] = combs #insert limited object
    combs[,nC+1] = signif(probVec2,signif) #round off probs for table:
    
    #Add to tables:
    table1 <- rbind(table1,cbind(Locus=loc,Rank=1:nrow(combs),combs))
    table1 <- rbind(table1, rep("",ncol(table1)) ) #add empty line
    table2 <- rbind(table2, c(toplist[[loc]]) )
    space <- cbind("","","","")
    for(cc in contrcn) {
      tmp <- deconvlistic[[loc]][[cc]]
      maxind <- maxI(p=tmp$Probability)
      newrow <- tmp[1:maxind,,drop=F]
      newrow[,2] <- signif(newrow[,2],signif)
      newrows <- as.matrix(cbind(cc,loc,newrow))
      table3 <- rbind(table3,newrows,space)
      
      tmp <- deconvlistica[[loc]][[cc]]
      maxind <- maxI2(p=tmp$Probability)
      newrow <- tmp[1:maxind,,drop=F]
      newrow[,2] <- signif(newrow[,2],signif)
      newrows <- as.matrix(cbind(cc,loc,newrow))
      table4 <- rbind(table4,newrows,space)
    }
    colnames(table3)[1:2] <- colnames(table4)[1:2] <- c("Contr.","Locus")
  } #end for each loci
  colnames(table2) <- paste0(topRankcn,"_",c(t(replicate(length( topRankcn),contrcn))))
  rownames(table2) <- locs
  
  #new version adds rankGic and rankGica (genotype ranks per contributor in addition to per allele)
  return(list(table1=table1,table2=table2,table3=table3,table4=table4,toprankGi=toplist,rankGi=deconvlisti))
} #end function

