#' @title deconvolve
#' @author Oyvind Bleka
#' @description deconvolve ranks the set of the most conditional posterior probability of genotypes the STR DNA mixture given a fitted model under a hypothesis.
#' @details The procedure calculates the likelihood for each single locus. Then it combines the most probable genotypes from each loci to produce a ranked list of deconvolved profiles.
#' 
#' @param mlefit Fitted object using calcMLE/contLikMLE function.
#' @param alpha Required sum of the listed posterior probabilities
#' @param maxlist NOT IN USE
#' @param signif Number of significant numbers
#' @param checkCalcs NOT IN USE
#' @return ret A list(table1,table2,table3,table4,rankGi,rankG,pG) where rankG is the ranked genotypes with corresponding probabilities in pG. rankgGi is the same, but per marker. table1 is now empty (since joint results is not obtained). table2 uses rankGi to find marginal results for top-genotypes. table3 and table4 shows this marginalized on genotypes and alleles per contributor 
#' @export

#signif=4;alpha=0.99
deconvolve = function(mlefit,alpha=0.99,maxlist=1000,signif=4,checkCalcs=FALSE){ 
  
  #Obtain list of marginal probs (already performed in MLE step)
  DCmargList = mlefit$DCmargList #get
  #if(is.null(DCmargList)) return(NULL)
  c <- mlefit$prepareC #returned from prepareC
  locs = c$markerNames #names(DCmargList)
  nC = c$nC #number of contrs

  #Helpfunctions to obtain marginal probabilities
  getMarg <- function(vals) { #get marginal of genotypes
    genos = names(vals)
    agg <- aggregate(vals,by=list(genos),sum) #get probabilities
    ord <- order(agg[,2],decreasing=TRUE)
    agg2 <- agg[ord,,drop=F]
    colnames(agg2) <- c("Genotype","Probability")
    return(agg2)
  }
  #Helpfunction to get marginal of alleles
  getMarg2 <- function(vals) { 
    genos = names(vals)
    tmp <- unlist(strsplit(genos,"/"))
    unA <- unique(tmp) #unique alleles
    x2 <- t(matrix(tmp,nrow=2))
    prob <- rep(NA,length(unA))  
    for(aa in unA) prob[which(unA==aa)] = sum(vals[rowSums(x2==aa)>0]) #sum probabilities
    ord <- order(prob,decreasing=TRUE)
    agg <- data.frame(Allele=unA[ord],Probability=prob[ord])
    return(agg)
  }

  #Step 2) Obtaining Deconvolved Genotype  results
  contrcn <- paste0("C",1:nC) #column name of contributors
  topRankcn <-  c("TopGenotype","probability","ratioToNextGenotype") #names for each contributor
  nTriAlleles = length(c$triAlleles)/3 #used to include tri-allele for knowns
  
  #OUTPUT VARIABLES
  #deconvlisti <- list() #stored full list (truncated on maxlist)
  deconvlistic <- list() #genotype list per contributor
  deconvlistica <- list() #allele list per contributor
  toplist <- list()
  table2 <- table3 <- table4 <- numeric()
  for(markerIdx in seq_along(c$markerNames)) { #extract info in c relevant for each markers:
# markerIdx=7;   print(markerIdx)
    loc = c$markerNames[markerIdx] #obtain locus name
    DCmarg <- DCmargList[[loc]] #obtain marginals
    Gset1 = paste0(c$genoList[[loc]][,1],"/",c$genoList[[loc]][,2])
    #Gset1 = colnames(DCmarg) #get genotypes
    nGenos = length(Gset1) #number of genos
    
    #obtain genotype of known individuals
    NOK0 = c$NOK 
    uind = 1:nC
    kind = integer()
    if(NOK0>0) {
      knownGind = c$knownGind[(markerIdx-1)*NOK0 + seq_len(NOK0)] #obtain genotype index for known contributors
      kind = which(knownGind != (-1)) #index of knowns (different from)
      uind = setdiff(uind,kind) #update unknown indice
    } 
    
    #Construct full matrix (also including knowns)
    DCmargFull = matrix(0,ncol=nC,nrow=nGenos,dimnames = list(Gset1,contrcn)) #construct full matrix
    if(length(uind)>0) DCmargFull[,uind] = t(DCmarg) #insert
    if(length(kind)>0) DCmargFull[cbind(knownGind[kind]+1,kind)] = 1 #insert certain genotype
    
     #Part 2: CALCULATE MARGINAL PRODUCTS (by-products of joint result deconvOBJ)
    #A) Calculate marginal probabilities for all contributors (genotypes and alleles):
    deconvlistica[[loc]] <- deconvlistic[[loc]] <- list()
    tab <- matrix(,nrow=3,ncol=nC)
    rownames(tab) <- topRankcn 
    colnames(tab) <- contrcn
    for(k in 1:nC) {
      DCmargVec = setNames(DCmargFull[,k],rownames(DCmargFull)) #this is more robust
      deconvlistic[[loc]][[contrcn[k]]] <- tmp <- getMarg(DCmargVec)
      deconvlistica[[loc]][[contrcn[k]]] <- getMarg2(DCmargVec)
      RTN <- ifelse(nrow(tmp)>1, signif(tmp[1,2]/tmp[2,2],signif),NA) #get ratio from first to second genotype
      if(is.infinite(RTN)) RTN = NA #set as NA
      tab[,k] <- c(tmp[1,1],signif(tmp[1,2],signif),RTN) 
    }
    toplist[[loc]] <- tab
    
    #Adding possible tri-alleles to table:
    if(nTriAlleles>0) {
      vectorElemsStartInds = 3*(seq_len(nTriAlleles)-1) + 1 #obtain element range in triAllele vector 
      markerInds = c$triAlleles[vectorElemsStartInds] + 1 #Obtain allele indices (long allelename vector)
      isThisMarker = markerInds==markerIdx
      if(any(isThisMarker)) { #if at least one element was within the range
        for(elemIdx in which(isThisMarker)) { #we traverse the elements for this marker
#          elemIdx = which(isThisMarker)[1]
          startIndAlleles = c$startIndMarker_nAlleles[markerIdx]
          triAlleleIdx =  c$triAlleles[vectorElemsStartInds[elemIdx]+1] + 1
          contrIdx =  c$triAlleles[vectorElemsStartInds[elemIdx]+2] + 1
          triAllele = c$alleleNames[startIndAlleles+triAlleleIdx]
          genotypeWithTriAllele = paste0(toplist[[loc]][1,contrIdx],"/",triAllele) #obtain added genotype
          toplist[[loc]][1,contrIdx] = genotypeWithTriAllele #UPDATING
          deconvlistic[[loc]][[contrcn[contrIdx]]][1,1] = genotypeWithTriAllele #UPDATING
          deconvlistica[[loc]][[contrcn[contrIdx]]] = rbind(c(triAllele,1),deconvlistica[[loc]][[contrcn[contrIdx]]]) #UPDATING
        }
      } 
    }
    
    #Add to tables:
    table2 <- rbind(table2, c(toplist[[loc]]) )
    space <- cbind("","","","")
    for(cc in contrcn) {
      #Obtain genotypes
      tmp <- deconvlistic[[loc]][[cc]]
      toprange <- which(cumsum(tmp$Probability)<=alpha)
      if(length(toprange)==0) toprange = 1
      newrow <- tmp[toprange,,drop=F]
      newrow[,2] <- signif(newrow[,2],signif)
      newrows <- as.matrix(cbind(cc,loc,newrow))
      table3 <- rbind(table3,newrows,space)
      
      #Obtain alleles
      tmp <- deconvlistica[[loc]][[cc]]
      toprange <-  which( tmp$Probability > (1-alpha)) #get most likely ones
      newrow <- tmp[toprange,,drop=F]
      newrow[,2] <- signif(as.numeric(newrow[,2]),signif)
      newrows <- as.matrix(cbind(cc,loc,newrow))
      table4 <- rbind(table4,newrows,space)
    }
    colnames(table3)[1:2] <- colnames(table4)[1:2] <- c("Contr.","Locus")
  } #end for each marker
  colnames(table2) <- paste0(topRankcn,"_",c(t(replicate(length( topRankcn),contrcn))))
  rownames(table2) <- c$markerNames
  
  #new version adds rankGic and rankGica (genotype ranks per contributor in addition to per allele)
  return(list(table1=NULL,table2=table2,table3=table3,table4=table4,toprankGi=toplist,rankGi=NULL))
} #end function

