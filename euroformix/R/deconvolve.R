#' @title deconvolve
#' @author Oyvind Bleka
#' @description deconvolve ranks the set of the most conditional posterior probability of genotypes the STR DNA mixture given a fitted model under a hypothesis.
#' @details The procedure calculates the likelihood for each single locus. Then it combines the most probable genotypes from each loci to produce a ranked list of deconvolved profiles.
#' 
#' Function calls procedure in c++ by using the package Armadillo and Boost.
#'
#' @param mlefit Fitted object using contLikMLE function.
#' @param alpha Required sum of the listed posterior probabilities.
#' @param maxlist The ranked deconvolved profile list will not exceed this number (used to avoid endless search).
#' @return ret A list(table1,table2,table3,table4,rankGi,rankG,pG) where rankG is the ranked genotypes with corresponding probabilities in pG. rankgGi is the same, but per marker. table1 is rankG and pG combined (joint results). table2 uses rankGi to find marginal results for top-genotypes. table3 and table4 shows this marginalized on genotypes and alleles per contributor 
#' @export
#' @references Cowell,R.G. et.al. (2014). Analysis of forensic DNA mixtures with artefacts. Applied Statistics, 64(1),1-32.
#' @keywords deconvolution
deconvolve = function(mlefit,alpha=0.95,maxlist=1000){
 theta2 <- theta <- mlefit$fit$thetahat #condition on mle parameter
 model <- mlefit$model #take out assumed model with given data
 locs <- names(model$popFreq)
 nL <- length(locs)
 nC <- model$nC
 xi <- model$xi
 fst <- model$fst
 np <- length(theta)#number of unknown parameters
 loglikYtheta <- function() {   #call c++- function: length(theta)=nC+1
   Cval  <- .C("loglikgammaC",as.numeric(0),as.numeric(theta),as.integer(np),ret$nC,ret$nK,ret$nL,ret$nS,ret$nA,ret$obsY,ret$obsA,ret$CnA,ret$allASind,ret$nAall,ret$CnAall,ret$Gvec,ret$nG,ret$CnG,ret$CnG2,ret$pG,ret$pA, as.numeric(model$prC), ret$condRef,as.numeric(model$threshT),fst,ret$mkvec,ret$nkval,as.numeric(model$lambda),ret$bp,as.integer(0),PACKAGE="euroformix")[[1]]
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
 dlist <- list()  #log-likelihood per genotype combinations per marker
 GClist2 <- list() #get index of ranked combined genotypes for all genotypes
 Glist2 <- list() #used to store names of genotypes 
 for(loc in locs) {
  samples <- lapply(model$samples,function(x) x[loc])
  ret <- prepareC(nC=nC,samples,popFreq=model$popFreq[loc],refData=model$refData[loc],condOrder=model$condOrder,knownRef=model$knownRef,kit=model$kit,model$knownRel,model$ibd,fst,incS=is.null(xi) || xi>0)
  uind <- which(ret$condRef==-1) #unknown genotype indices
  kind <- which(ret$condRef!=-1) #known genotype indices
  nU <- length(uind) #number of unknowns
  if(nU==0) { #if no unknowns
    dlist[[loc]] <- 0
    GClist2[[loc]] <- t(as.matrix(ret$condRef + 1)) #insert the conditioned ones
    Glist2[[loc]] <- calcGjoint(freq=model$popFreq[[loc]],nU=1)$G #Get genotypes
    next
  }
  #CALCULATING JOINT PROB GENOTYPE FOR ALL COMBINATIONS:
  #MUST SPECIFY MODEL OF model object here
  pGlist = calcGjoint(freq=model$popFreq[[loc]],nU=nU,fst = model$fst, refK = unlist(model$refData[[loc]][model$knownRef]), refR = unlist(model$refData[[loc]][model$knownRel]), ibd = model$ibd, sortComb = TRUE)
  Glist2[[loc]] =  pGlist$G #store original genotype outcomes

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
   genoUind = combGind[gind,] #get row
   ret$condRef[uind] <- as.integer(genoUind - 1) #insert genotypes to consider
   dvec[gind] <- loglikYtheta() #insert log P(E|G)

   if(nU>2) genoUind[-1] <- sort(genoUind[-1],decreasing=TRUE) #sort the indices to be ordered from 2.index
   pg = pGlist$Gprob[ rbind(genoUind) ]  #CALCULATING JOINT PROB GENOTYPE
   dvec[gind] <- dvec[gind] + log(pg) #insert log joint prob
 # dvec[gind] <- dvec[gind] + sum(log(ret$pG[combGind[gind,]])) #PREV VERSION
  }
  combGind2 <- numeric()
  for(c in 1:nC) {
   if(c%in%kind) combGind2 <- cbind(combGind2, rep(ret$condRef[c]+1,nrow(combGind)) ) #add genotype index of the reference
   if(c%in%uind) combGind2 <- cbind(combGind2, combGind[,which(c==uind)]) #add genotype index of the reference
  }
  isOK <- !is.infinite(dvec) #remove genotypes giving zero likelihood
  combGind2 <- combGind2[isOK,,drop=F]
  dvec <- dvec[isOK]
  rank <- order(dvec,decreasing=TRUE)
  dlist[[loc]] <- dvec[rank]   #log-likelihood per per marker per
  GClist2[[loc]] <- combGind2[rank,,drop=F] #get index of ranked combined genotypes 
 }

 kvec <- 1:nC
 colN <- paste0("C",kvec) #column name

 #Step 3) Convert rank-list to list with allele-names
 deconvlisti <- list() #list per locus for all contributors
 for(loc in locs) {
   ii <- which(locs==loc)
   genv <-  paste0(Glist2[[loc]][,1],"/",Glist2[[loc]][,2]) #get vector of genotypes
   pGi <- exp(dlist[[loc]]) #convert to normal scale
   deconvlisti[[loc]] <- matrix(genv[ GClist2[[loc]] ],nrow=nrow(GClist2[[loc]])) #translate indices to genotype names
   deconvlisti[[loc]] <- cbind(deconvlisti[[loc]], pGi/sum(pGi) ) #add probabilities per makers
   colnames(deconvlisti[[loc]]) <- c(colN,"Probability")
 }

 #Step4) Create table layouts:

 #A) Calculate marginal probabilities for all contributors (genotypes and alleles):
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
 deconvlistic <- list() #genotype list per contributor
 deconvlistica <- list() #allele list per contributor
 cn <-  c("TopGenotype","probability","ratioToNextGenotype") #names for each contributor
 toplist <- list()
 for(loc in locs) {
  deconvlistica[[loc]] <- deconvlistic[[loc]] <- list()
  X <- deconvlisti[[loc]]
  nc <- ncol(X) #number of column
  nr <- nc - 1 #number of contributors
  tab <- matrix(,nrow=3,ncol=nr)
  rownames(tab) <- cn 
  colnames(tab) <- colN
  for(rr in 1:nr) {
   deconvlistic[[loc]][[colN[rr]]] <- tmp <- getMarg(x=X[,rr],y=as.numeric(X[,nc]))
   deconvlistica[[loc]][[colN[rr]]] <- getMarg2(x=X[,rr],y=as.numeric(X[,nc]))
   rat <- ifelse(nrow(tmp)>1, tmp[1,2]/tmp[2,2],NA) #get ratio from first to second genotype
   tab[,rr] <- c(tmp[1,1],signif(tmp[1,2],4),signif(rat,4)) 
  }
  toplist[[loc]] <- tab
 }
 
 maxI <- function(p) min(min(which(cumsum(p)>=alpha)),maxlist,length(p))
 #B) Create tables
 table1 <- table2 <- table3 <- table4 <- numeric()
 for(loc in locs) {
   combs <- deconvlisti[[loc]]
   prob <- as.numeric(combs[,ncol(combs)])
   combs[,ncol(combs)] <- signif(prob,4)
   maxind <- maxI(prob)
   combs <- combs[1:maxind,,drop=F]
   table1 <- rbind(table1,cbind(loc,1:nrow(combs),combs))
   table1 <- rbind(table1, rep("",ncol(table1)) )
 }
 colnames(table1)[1:2] <- c("Locus","Rank")

 for(loc in locs) table2 <- rbind(table2, c(toplist[[loc]]) )
 colnames(table2) <- paste0(cn,"_",c(t(replicate(length(cn),colN))))
 rownames(table2) <- locs

 maxI2 <- function(p) min(max(which(p>(1-alpha))),maxlist,length(p))
 space <- cbind("","","","")
 for(cc in colN) {
  for(loc in locs) {
   tmp <- deconvlistic[[loc]][[cc]]
   maxind <- maxI(p=tmp$Probability)
   newrow <- tmp[1:maxind,,drop=F]
   newrow[,2] <- signif(newrow[,2],4)
   newrows <- as.matrix(cbind(cc,loc,newrow))
   table3 <- rbind(table3,newrows,space)

   tmp <- deconvlistica[[loc]][[cc]]
   maxind <- maxI2(p=tmp$Probability)
   newrow <- tmp[1:maxind,,drop=F]
   newrow[,2] <- signif(newrow[,2],4)
   newrows <- as.matrix(cbind(cc,loc,newrow))
   table4 <- rbind(table4,newrows,space)
  }
 }
 colnames(table3)[1:2] <- colnames(table4)[1:2] <- c("Contr.","Locus")

 #new version adds rankGic and rankGica (genotype ranks per contributor in addition to per allele)
 return(list(table1=table1,table2=table2,table3=table3,table4=table4,toprankGi=toplist,rankGi=deconvlisti))
#, rankGic=deconvlistic,rankGica=deconvlistica)
} #end function

