#' @title deconvolve
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description deconvolve ranks the set of the most conditional posterior probability of genotypes the STR DNA mixture given a fitted model under a hypothesis.
#' @details The procedure calculates the likelihood for each single locus. Then it combines the most probable genotypes from each loci to produce a ranked list of deconvolved profiles.
#' 
#' Function calls procedure in c++ by using the package Armadillo and Boost.
#'
#' @param mlefit Fitted object using contLikMLE function.
#' @param alpha Required sum of the listed posterior probabilities.
#' @param maxlist The ranked deconvolved profile list will not exceed this number (used to avoid endless search).
#' @param unknownonly A boolean whether table should only contain the unknown genotypes or both known and unknown genotypes.
#' @return ret A list(table1,table2,rankGi,rankG,pG) where rankG is the ranked genotypes with corresponding probabilities in pG. rankgGi is the same, but per marker. table1 is rankG and pG combined (joint results). table2 uses rankGi to find marginal results for top-genotypes.
#' @export
#' @references Cowell,R.G. et.al. (2014). Analysis of forensic DNA mixtures with artefacts. Applied Statistics, 64(1),1-32.
#' @keywords deconvolution
deconvolve = function(mlefit,alpha=0.95,maxlist=1000,unknownonly=TRUE){
 theta2 <- theta <- mlefit$fit$thetahat #condition on mle parameter
 model <- mlefit$model #take out assumed model with given data
 locs <- names(model$popFreq)
 nL <- length(locs)
 nC <- model$nC
 xi <- model$xi
 np <- length(theta)#number of unknown parameters
 loglikYtheta <- function() {   #call c++- function: length(theta)=nC+1
   Cval  <- .C("loglikgammaC",as.numeric(0),as.numeric(theta),as.integer(np),ret$nC,ret$nK,ret$nL,ret$nS,ret$nA,ret$obsY,ret$obsA,ret$CnA,ret$allAbpind,ret$nAall,ret$CnAall,ret$Gvec,ret$nG,ret$CnG,ret$CnG2,ret$pG,ret$pA, as.numeric(model$prC), ret$condRef,as.numeric(model$threshT),as.numeric(model$fst),ret$mkvec,ret$nkval,as.numeric(model$lambda),ret$bp,as.integer(0),PACKAGE="euroformix")[[1]]
   return(Cval) #weight with prior of tau and 
 }
 nodeg <- is.null(model$kit) #check for degradation
 if(nodeg) theta <- c(theta[1:(nC+1)],1) #insert beta variable equal 1
 if(!is.null(xi)) {
  theta <- c(theta,as.numeric(xi))
 } else {
  if(nodeg) theta <- c(theta,theta2[np])
 }
 #Using information in ret to try out different genotypes:

 #Step 1) Calculate L(E|g,thetahat) for each marker
 dlist <- list()
 GClist <- list()
 for(loc in locs) {
  samples <- lapply(model$samples,function(x) x[loc])
  ret <- prepareC(nC=nC,samples,popFreq=model$popFreq[loc],refData=model$refData[loc],condOrder=model$condOrder,knownRef=model$knownRef,kit=model$kit)
  uind <- which(ret$condRef==-1) #unknown genotype indices
  nU <- length(uind) #number of unknowns
  if(nU==0) stop("There was no unknown genotype profiles to estimate. The evaluation will not be done!")
  ret$nK <- nC  #number of known is equal number of contributors
  Gset <- matrix(ret$Gvec,ncol=2) #genotype possibilities
  Glist <- list()
  for(k in 1:nU) {
   Glist[[k]] <- 1:nrow(Gset)
  }
  combGind <- expand.grid(Glist) #get all combinations
  combGind <- as.matrix(combGind,nrow=nrow(combGind))
  dvec <- rep(NA,nrow(combGind))
  #calculate for each genotypes:
  for(gind in 1:nrow(combGind)) { #for each possible genotypes:
   ret$condRef[uind] <- as.integer(combGind[gind,] - 1) #genotypes to consider
   dvec[gind] <- loglikYtheta() 
   dvec[gind] <- dvec[gind] + sum(log(ret$pG[combGind[gind,]])) #add genotype probability as well
  }
  isOK <- !is.infinite(dvec)
  combGind <- combGind[isOK,]
  dvec <- dvec[isOK]
  rank <- order(dvec,decreasing=TRUE)
  dlist[[loc]] <- dvec[rank] 

  if(is.null(dim(combGind))) { #threat the case of one unknown
   if(nU==1) GClist[[loc]] <- matrix(combGind[rank],ncol=1) #it's a vector because of only 1 unknown
   if(nU>1) GClist[[loc]] <- matrix(combGind,nrow=1) #it's a vector since it was only 1 possibility
  } else {
   GClist[[loc]] <- combGind[rank,]
  }
 }

 #Step 2) Combine markers to create full profiles (this is it's own function):
 loghdval <- mlefit$fit$loglik 
 #loghdvali <- logLiki(mlefit) #Not necessary: get likelihood for each marker: 
 if(is.null(xi)) loghdval <- loghdval - log(model$pXi(theta2[np])) #reported likval has taken pXi into account
 rankGlist <- combineRank(dlist,loghdval=loghdval,alpha=alpha,maxsearch=maxlist)
 pG <- rankGlist$pG
 Gset  <- rankGlist$rankG

 kvec <- 1:nC
 if(unknownonly) kvec <- uind  

 #Step 3) Convert rank-list to list with allele-names
 Glist <- getGlist(model$popFreq) #get genotype list with genotypes and corresponding frequencies
 deconvlist <- list()
 deconvlisti <- list() #list per locus
 for(i in 1:nL) { #convert names for each locus
  rankgeno <- GClist[[locs[i]]][Gset[,i],]
  if(is.null(dim(rankgeno))) rankgeno <- as.matrix(rankgeno) #make matrix again
  rankgenos <- numeric()
  rankgenos2<- numeric() #added: per locus
  for(k in 1:nC) { #for each contributor
   if(k%in%uind) { #if unknown contributor
    geno <- matrix(Glist[[locs[i]]]$G[rankgeno[,which(k==uind)],],ncol=2) #get allele named genotype
    geno2 <- matrix(Glist[[locs[i]]]$G[GClist[[locs[i]]][,which(k==uind)],],ncol=2) #get allele named genotype
   } else if(!unknownonly) { #if known contributors in addition(they are given in reference) 
    geno <- sort(model$refData[[locs[i]]][[which(model$condOrder==k)]])
    geno <- matrix( rep(geno,nrow(rankgeno)),ncol=2,byrow=TRUE)

    geno2 <- sort(model$refData[[locs[i]]][[which(model$condOrder==k)]])
    geno2 <- matrix( rep(geno2,nrow(GClist[[locs[i]]])),ncol=2,byrow=TRUE)
   } else {
    next #skip to next contributor
   }
   geno <- paste0(geno[,1],"/",geno[,2])
   rankgenos <- cbind(rankgenos,geno)

   geno2 <- paste0(geno2[,1],"/",geno2[,2])
   rankgenos2 <- cbind(rankgenos2,geno2)
  }
  colnames(rankgenos) <- colnames(rankgenos2) <- paste0("C",kvec)
  deconvlist[[locs[i]]] <- rankgenos 

  #ADDED: create per-locus probabilities. Enough to normalize over all possibilities
  pGi <- exp(dlist[[locs[i]]])
  deconvlisti[[locs[i]]] <- cbind(rankgenos2, pGi/sum(pGi))#exp(dlist[[locs[i]]]-loghdvali[i]))
 }

 #Step4) Create table layouts:

 #A) Create top-genotype probabilities for all unknown contr:
 cn <-  c("TopGenotype","probability","ratioToNextGenotype") #names for each unknown
 getMarg <- function(x,y) {
  agg <- aggregate(y,by=list(x),sum) #get probabilities
  ord <- order(agg[,2],decreasing=TRUE)
  agg2 <- agg[ord,]
  return(agg2)
 }
 toplist <- lapply(deconvlisti,function(X) {
  nc <- ncol(X) #number of column
  nr <- nc - 1 #number of unknown
  tab <- matrix(,nrow=3,ncol=nr)
  rownames(tab) <- cn 
  colnames(tab) <- paste0("unknown",1:nr)
  for(rr in 1:nr) {
   agg <- getMarg(x=X[,rr],y=as.numeric(X[,nc]))
   rat <- NA
   if(nrow(agg)>1) rat <- agg[1,2]/agg[2,2] #get ratio from first to second 
   print(agg)
   tab[,rr] <- c(agg$Group.1[1],signif(agg$x[1],4),signif(rat,4)) 
  }
  return(tab)
 })

 #B) Create tables
 table1 <- table2 <- numeric()
 for(loc in locs) table1 <- cbind(table1,deconvlist[[loc]])
 table1 <-  cbind(table1,signif(pG,4))
 colnames(table1) <- c(paste0(c(t(replicate(length(kvec),locs))),"_C",kvec),"probability")

 for(loc in locs) table2 <- rbind(table2, c(toplist[[loc]]) )
 colnames(table2) <- paste0(cn,"_C",c(t(replicate(length(cn),kvec))))
 rownames(table2) <- locs

 return(list(table1=table1,table2=table2,toprankGi=toplist,rankG=deconvlist,rankGi=deconvlisti,pG=pG))
} #end function
