
#' @title contLikMLEpara
#' @author Oyvind Bleka
#' @description Same as the contLikMLE function (dummy function): earlier version did parallelization
#'
#' @param nC Number of contributors in model. Must be a constant.
#' @param samples A List with samples which for each samples has locus-list elements with list elements adata and hdata. 'adata' is a qualitative (allele) data vector and 'hdata' is a quantitative (peak heights) data vector.
#' @param popFreq A list of allele frequencies for a given population.
#' @param refData Reference objects with (2-size) allele-vector given in list element [[i]][[s]].
#' @param condOrder Specify conditioning references from refData (must be consistent order). For instance condOrder=(0,2,1,0) means that we restrict the model such that Ref2 and Ref3 are respectively conditioned as 2. contributor and 1. contributor in the model. 
#' @param knownRef Specify known non-contributing references from refData (index). For instance knownRef=(1,2) means that reference 1 and 2 is known non-contributor in the hypothesis. This affectes coancestry correction.
#' @param xi A numeric giving stutter-ratio if it is known. Default is 0, meaning stutter is not used.
#' @param prC A numeric for allele drop-in probability. Default is 0.
#' @param nDone Number of optimizations required providing equivalent results (same logLik value obtained)
#' @param threshT The detection threshold given. Used when considering probability of allele drop-outs.
#' @param fst The co-ancestry coefficient. Default is 0.
#' @param lambda Parameter in modeled peak height shifted exponential model. Default is 0.
#' @param pXi Prior function for xi-parameter (stutter). Flat prior on [0,1] is default.
#' @param delta Scaling of variation of normal distribution when drawing random startpoints. Default is 1.
#' @param kit Used to model degradation. Must be one of the shortnames of kit: {"ESX17","ESI17","ESI17Fast","ESX17Fast","Y23","Identifiler","NGM","ESSPlex","ESSplexSE","NGMSElect","SGMPlus","ESX16", "Fusion","GlobalFiler"}. 
#' @param verbose Whether printing optimization progress. Default is TRUE.
#' @param maxIter Maximum number of iterations for the optimization. Default is 100.
#' @param knownRel gives the index of the reference which the 1st unknown is related to.
#' @param ibd the identical by decent coefficients of the relationship (specifies the type of relationship)
#' @param xiFW A numeric giving FW stutter-ratio if it is known.Default is 0, meaning stutter is not used.
#' @param pXiFW Prior function for xiFW-parameter (FW stutter). Flat prior on [0,1] is default.
#' @param seed The user can set seed if wanted
#' @param maxThreads Maximum number of threads to be executed by the parallelization
#' @param steptol Argument used in the nlm function for faster return from the optimization (tradeoff is lower accuracy).
#' @return ret A list(fit,model,nDone,delta,prepareC) where fit is Maximixed likelihood elements for given model.
#' @export

contLikMLEpara = function(nC,samples,popFreq,refData=NULL,condOrder=NULL,knownRef=NULL,xi=0,prC=0,nDone=2,threshT=50,fst=0,lambda=0,pXi=function(x)1,delta=1,kit=NULL,verbose=TRUE,maxIter=100,knownRel=NULL,ibd=c(1,0,0),xiFW=0,pXiFW=function(x)1,seed=NULL,maxThreads=32,steptol=1e-3){

 #call function contLikMLE with same arguments
 mle = euroformix::contLikMLE(nC=nC,samples=samples,popFreq,refData=refData,condOrder=condOrder,knownRef=knownRef,xi=xi,prC=prC,nDone=nDone,threshT=threshT,fst=fst,lambda=lambda,pXi=pXi,delta=delta,kit=kit,verbose=verbose,maxIter=maxIter,knownRel=knownRel,ibd=ibd,xiFW=xiFW,pXiFW=pXiFW,seed=seed,maxThreads=maxThreads,steptol=steptol)
 return(mle)

} #end function

