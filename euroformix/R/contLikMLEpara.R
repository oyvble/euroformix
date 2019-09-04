
#' @title contLikMLEpara
#' @author Oyvind Bleka
#' @description Parallelization on contLikMLE using snow
#' @details The procedure is doing parallelization of the contLikMLE function
#' 
#' The procedure also does a Laplace Approximation of the marginalized likelihood (theta integrated out) and returns the log-marginal likelihood as logmargL in the fit-list.
#' 
#' Function calls procedure in c++ by using the package Armadillo and Boost.
#'
#' @param nC Number of contributors in model.
#' @param samples A List with samples which for each samples has locus-list elements with list elements adata and hdata. 'adata' is a qualitative (allele) data vector and 'hdata' is a quantitative (peak heights) data vector.
#' @param popFreq A list of allele frequencies for a given population.
#' @param refData Reference objects with (2-size) allele-vector given in list element [[i]][[s]].
#' @param condOrder Specify conditioning references from refData (must be consistent order). For instance condOrder=(0,2,1,0) means that we restrict the model such that Ref2 and Ref3 are respectively conditioned as 2. contributor and 1. contributor in the model. 
#' @param knownRef Specify known non-contributing references from refData (index). For instance knownRef=(1,2) means that reference 1 and 2 is known non-contributor in the hypothesis. This affectes coancestry correction.
#' @param xi A numeric giving stutter-ratio if it is known. Default is NULL, meaning it is integrated out.
#' @param prC A numeric for allele drop-in probability. Default is 0.
#' @param nDone Maximum number of random evaluations nlm-optimizing routing. Default is 1.
#' @param threshT The detection threshold given. Used when considering probability of allele drop-outs.
#' @param fst The co-ancestry coefficient. Default is 0.
#' @param lambda Parameter in modeled peak height shifted exponential model. Default is 0.
#' @param pXi Prior function for xi-parameter (stutter). Flat prior on [0,1] is default.
#' @param delta Standard deviation of normal distribution when drawing random startpoints. Default is 10.
#' @param kit Used to model degradation. Must be one of the shortnames of kit: {"ESX17","ESI17","ESI17Fast","ESX17Fast","Y23","Identifiler","NGM","ESSPlex","ESSplexSE","NGMSElect","SGMPlus","ESX16", "Fusion","GlobalFiler"}. 
#' @param verbose Boolean whether printing optimization progress. Default is TRUE.
#' @param maxIter Maximum number of iterations for the optimization. Default is 30.
#' @param knownRel gives the index of the reference which the 1st unknown is related to.
#' @param ibd the identical by decent coefficients of the relationship (specifies the type of relationship)
#' @return ret A list(fit,model,nDone,delta) where fit is Maximixed likelihood elements for given model.
#' @export
#' @references Cowell,R.G. et.al. (2014). Analysis of forensic DNA mixtures with artefacts. Applied Statistics, 64(1),1-32.
#' @keywords continuous model, Maximum Likelihood Estimation

contLikMLEpara = function(nC,samples,popFreq,refData=NULL,condOrder=NULL,knownRef=NULL,xi=NULL,prC=0,nDone=1,threshT=50,fst=0,lambda=0,pXi=function(x)1,delta=10,kit=NULL,verbose=TRUE,maxIter=100,knownRel=NULL,ibd=c(1,0,0)){
 library(parallel, warn.conflicts = FALSE)

 ncores <- parallel::detectCores() #number of physical cores (parallel)
 nCl <- min(ncores,nDone) #number of clusters
 print(paste0("\nNumber of optimalizations will be ",nCl ))
 #if( nCl==1 || nCl%%2!=0) stop("Please change the number of startpoints to an even number. This is can be changed under Optimization.")
 cl <- parallel::makeCluster(nCl,type="PSOCK")
 inputlist <- list(nC=nC,samples=samples,popFreq=popFreq,refData=refData,condOrder=condOrder,knownRef=knownRef,xi=xi,prC=prC,nDone=1,threshT=threshT,fst=fst,lambda=lambda,pXi=pXi,delta=delta,kit=kit,verbose=verbose,maxIter=maxIter,knownRel=knownRel,ibd=ibd)
 inputlist <- rep(list(inputlist),nCl) #number of clusters

 fclust<-function(x) {
  library(euroformix)
  return(euroformix::contLikMLE(x$nC,x$samples,x$popFreq,x$refData,x$condOrder,x$knownRef,x$xi,x$prC,x$nDone,x$threshT,x$fst,x$lambda,x$pXi,x$delta,x$kit,x$verbose,x$maxIter,x$knownRel,x$ibd))
 }
 retlist <- parallel::clusterApply(cl,inputlist,fclust)
 parallel::stopCluster(cl)
 ismax <- which.max(sapply(retlist,function(x) x$fit$loglik))[1]
 return(retlist[[ismax]])
} #end function

