#' @title deconvolve
#' @author Oyvind Bleka
#' @description deconvolve ranks the set of the most conditional posterior probability of genotypes the STR DNA mixture given a fitted model under a hypothesis.
#' @details The procedure calculates the likelihood for each single locus. Then it combines the most probable genotypes from each loci to produce a ranked list of deconvolved profiles.
#' 
#' Function calls the likelihood procedure in C++ by using the package Boost.
#'
#' @param mlefit Fitted object using contLikMLE function.
#' @param alpha Required sum of the listed posterior probabilities.
#' @param maxlist The ranked deconvolved profile list will not exceed this number (used to avoid endless search).
#' @param verbose Boolean whether printing deconvolution progress. Default is TRUE.
#' @return ret A list(table1,table2,table3,table4,rankGi,rankG,pG) where rankG is the ranked genotypes with corresponding probabilities in pG. rankgGi is the same, but per marker. table1 is rankG and pG combined (joint results). table2 uses rankGi to find marginal results for top-genotypes. table3 and table4 shows this marginalized on genotypes and alleles per contributor 
#' @export
#' @references Cowell,R.G. et.al. (2014). Analysis of forensic DNA mixtures with artefacts. Applied Statistics, 64(1),1-32.
#' @keywords deconvolution
deconvolve = function(mlefit,alpha=0.95,maxlist=1000,verbose=TRUE){ #alpha=0.95;maxlist=1000;verbose=TRUE
 thhat <- mlefit$fit$thetahat #condition on mle parameter
 model <- mlefit$model #take out assumed model with given data

 c <- mlefit$prepareC #returned from prepareC
 locs = c$locs #will be same as names(model$popFreq) #get loci names to evaluate
 nC = c$nC #number of contributors
 usedeg <- !is.null(model$kit) #check for degradation model
 
 #Prepare fixed params:
 nM = c$nM #number of markers to evaluate
 ATv = model$threshT
 pCv = model$prC
 lambdav = model$lambda
 fstv = model$fst
 
 #Prepare estimated params:
 mixprop = head(thhat,nC-1)# fitted mix proportions
 mu1 = thhat[nC] #fitted PHexp
 sig1 = thhat[nC+1] #fitted PHvar
 
 #Obtain remaining params
 beta1 = 1 #set default values
 xiB = model$xi #set default values
 xiF = model$xiFW #set default values
 indsel = nC + 2  #get index to insert estimated variable form param vector(thhat)
 if(!is.null(model$kit)) {
   beta1 = thhat[indsel] #get degrad param
   indsel = indsel + 1 #update index
 }
 if(is.null(xiB)) {
   xiB = thhat[indsel] #get BWstutter param
   indsel = indsel + 1 #update index
 }
 if(is.null(xiF)) xiF = thhat[indsel] #get FWstutter param
 
 #Step 1) Calculate L(E|g,thetahat) for each marker
 dlist <- list()  #joint genotype probs for each combinations per marker
 GClist2 <- list() #get index of ranked combined genotypes for all genotypes
 Glist2 <- list() #used to store names of genotypes 
 
 startIndPS = 0; #start marker index for potential stutters (own vectors)
 startIndMarker1=0;#start marker index (1 rep)
 startIndMarker2=0;#start marker index (nRep[m] reps)
 
 progcount  = 1 #progress counter
 if(verbose) progbar <- txtProgressBar(min = 0, max =nM, style = 3) #create progress bar
 
 for(m in 1:nM) { #extract info in c relevant for each markers:
  loc = locs[m]
  ind1 = startIndMarker1 + 1:c$nA[m] #get index of 1 repitition
  ind2 = startIndMarker2 + 1:(c$nA[m]*c$nRep[m]) #get index of 1 repitition
  indPS = startIndPS +  1:c$nPS[m] #get index of potential stutters
  
  #Update indices for next marker (must be done early):
  startIndMarker1 = startIndMarker1 + c$nA[m] #get start position of marker m+1 (1 rep)
  startIndMarker2 = startIndMarker2 + c$nA[m]*c$nRep[m]; #get start position of marker m+1 (nRep), used only for Peaks only
  startIndPS = startIndPS + c$nPS[m]; #add number of potential stutters
  
  if(c$nPS[m]==0) indPS = numeric()
  indKnownGind = (nC*(m-1)+1):(nC*m) #get index where the genotype indices are
  knownGind = c$knownGind[indKnownGind] #extract the vector with genotype indices (-1 means unknown)
  
  uind <- which(knownGind==-1) #unknown genotype indices
  kind <- which(knownGind!=-1) #known genotype indices
  nU <- length(uind) #number of unknowns
  if(nU==0) { #if no unknowns
    dlist[[loc]] <- 0
    GClist2[[loc]] <- t(as.matrix(knownGind + 1)) #insert the conditional genotypes
    Glist2[[loc]] <- calcGjoint(freq=model$popFreq[[loc]],nU=1)$G #Get genotypes
    next #skip to next marker
  }
  #CALCULATING JOINT PROB GENOTYPE FOR ALL COMBINATIONS:
  #MUST SPECIFY MODEL OF model object here
  if(!all(model$popFreq[[loc]]==c$FvecLong[ind1])) stop("Frequencies in model$popFreq and prepareC was different!") #extract frequencies
  refInd = c(which(model$condOrder>0),model$knownRef) #obtain ref index of BOTH conditional contr and known non-contributors (must be given as known ref alleles)
  pGlist = calcGjoint(freq=model$popFreq[[loc]],nU=nU,fst = fstv[m], refK = unlist(model$refData[[loc]][refInd]) , refR = unlist(model$refData[[loc]][model$knownRel]), ibd = model$ibd, sortComb = TRUE)
  Gset <- Glist2[[loc]] <- pGlist$G #genotype possibilities

  Glist <- list() #get genotype index list of all outcome
  for(k in 1:nU) {
    Glist[[k]] <- 1:nrow(Gset)
  }
  combGind <- expand.grid(Glist) #get all combinations
  combGind <- as.matrix(combGind,nrow=nrow(combGind))
  dvec <- rep(NA,nrow(combGind))
  
  #calculate for each genotypes:
  for(gind in 1:nrow(combGind)) { #for each possible genotypes:
    genoUind = combGind[gind,] #get row
    knownGind2 = knownGind 
    knownGind2[uind] =  as.integer(genoUind - 1) #insert genotypes to consider
#    Gset[knownGind2+1,] #show Genotype combinations
    dvec[gind]  = .C("loglikgammaC",as.numeric(0),as.integer(nC),as.integer(nC),as.integer(knownGind2),as.numeric(mixprop),as.numeric(mu1),as.numeric(sig1),as.numeric(beta1),as.numeric(xiB),as.numeric(xiF),as.numeric(ATv[m]),as.numeric(pCv[m]),as.numeric(lambdav[m]),as.numeric(fstv[m]),c$nReps[m],as.integer(1),c$nA[m],c$YvecLong[ind2],c$FvecLong[ind1],c$nTypedLong[m],c$maTypedLong[ind1],c$basepairLong[ind1],c$BWvecLong[ind1],c$FWvecLong[ind1],c$nPS[m],c$BWPvecLong[indPS],c$FWPvecLong[indPS],as.integer(1),as.integer(0),as.integer(0),c$relGind,c$ibdLong,PACKAGE="euroformix")[[1]] #Note: genotype probailities (related) is not considered here, since the prob is calculated after

    if(nU>2) genoUind[-1] <- sort(genoUind[-1],decreasing=TRUE) #sort the indices to be ordered from 2.index
    pg = pGlist$Gprob[ rbind(genoUind) ]  #CALCULATING JOINT PROB GENOTYPE (may include relatedness probs)
    dvec[gind] <- dvec[gind] + log(pg) #insert log joint prob
  }
  
  combGind2 <- numeric() #including the genotype of all contributors:
  for(k in 1:nC) { #for each contributors
    if(k%in%kind) combGind2 <- cbind(combGind2, rep(knownGind[k]+1,nrow(combGind)) ) #add genotype index of the reference
    if(k%in%uind) combGind2 <- cbind(combGind2, combGind[,which(k==uind)]) #add genotype index of the reference
  }
  isOK <- !is.infinite(dvec) #remove genotypes giving zero likelihood
  combGind2 <- combGind2[isOK,,drop=F] #removing genotypes which gives zero likelihood 
  dvec <- dvec[isOK]
  rank <- order(dvec,decreasing=TRUE) #order the genotypes wrt post prob values
  dlist[[loc]] <- dvec[rank]   #log-likelihood per per marker per
  GClist2[[loc]] <- combGind2[rank,,drop=F] #get index of ranked combined genotypes 
  
  if(verbose) { #only show progressbar if verbose (update)
    setTxtProgressBar(progbar,progcount)
    progcount <- progcount + 1
  } 
  
 } #end for each loci
 
 #POST PROCESSING:
 kvec <- 1:nC #index of contributors
 colN <- paste0("C",kvec) #column name of contributors

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
 maxI <- function(p) min(min(which(cumsum(p)>=alpha)),maxlist,length(p))  #helpfunction to obtain a maximum size of a vector (bounded in both length and probability)
 
 #A) Calculate marginal probabilities for all contributors (genotypes and alleles):
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

