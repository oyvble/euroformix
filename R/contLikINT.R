
#' @title contLikINT
#' @author Oyvind Bleka
#' @description contLikINT marginalizes the likelihood through numerical integration.
#' @details Replaced by new function calcINT
#' 
#' @param nC Number of contributors in model.
#' @param samples A List with samples which for each samples has locus-list elements with list elements adata and hdata. 'adata' is a qualitative (allele) data vector and 'hdata' is a quantitative (peak heights) data vector.
#' @param popFreq A list of allele frequencies for a given population.
#' @param lower Lower bounds of parameters. Must be in following order: mx1,..,mx_(nC-1),mu,sigma,beta,xi.
#' @param upper Upper bounds of parameters. Must be in following order: mx1,..,mx_(nC-1),mu,sigma,beta,xi.
#' @param refData Reference objects has locus-list element [[i]] with a list element 'r' which contains a 2 long vector with alleles for each references.
#' @param condOrder Specify conditioning references from refData (must be consistent order). For instance condOrder=(0,2,1,0) means that we restrict the model such that Ref2 and Ref3 are respectively conditioned as 2. contributor and 1. contributor in the model. 
#' @param knownRef Specify known non-contributing references from refData (index). For instance knownRef=(1,2) means that reference 1 and 2 is known non-contributor in the hypothesis. This affectes coancestry correction.
#' @param xi A numeric giving stutter-ratio if it is known. Default is NULL, meaning it is integrated out.
#' @param prC A numeric for allele drop-in probability. Default is 0.
#' @param reltol Required relative tolerance error of evaluations in integration routine. 
#' @param threshT The detection threshold given. Used when considering probability of allele drop-outs.
#' @param fst is the coancestry coeffecient. Default is 0.
#' @param lambda Parameter in modeled peak height shifted exponential model. Default is 0.
#' @param pXi Prior function for xi-parameter (stutter). Flat prior on [0,1] is default.
#' @param kit shortname of kit: Obtained from getKit()
#' @param scale used to make integrale calculateable for small numbers. For scale!=0, integrale must be scaled afterwards with exp(-scale) to be correct.
#' @param maxEval Maximum number of evaluations in the adaptIntegrate function. Default is 0 which gives an infinite limit.
#' @param knownRel gives the index of the reference which the 1st unknown is related to.
#' @param ibd the identical by decent coefficients of the relationship (specifies the type of relationship)
#' @param xiFW A numeric giving FW stutter-ratio if it is known.Default is 0, meaning stutter is not used.
#' @param pXiFW Prior function for xiFW-parameter (FW stutter). Flat prior on [0,1] is default.
#' @param maxThreads Maximum number of threads to be executed by the parallelization
#' @param verbose Whether printing limits to integrate over. Printing progress if maxEval>0. Default is TRUE.
#' @return ret A list(margL,deviation,nEvals,scale) where margL is Marginalized likelihood for hypothesis (model) given observed evidence, deviation is the confidence-interval of margL, nEvals is number of evaluations.
#' @export 
#' @references Hahn,T. (2005). CUBA - a library for multidimensional numerical integration. Computer Physics Communications, 168(2),78-95.
#' @keywords Marginalized likelihood


contLikINT = function(nC,samples,popFreq,lower,upper,refData=NULL,condOrder=NULL,knownRef=NULL,xi=NULL,prC=0,reltol=1,threshT=50,fst=0,lambda=0,pXi=function(x)1,kit=NULL,scale=0,maxEval=0,knownRel=NULL,ibd=c(1,0,0),xiFW=0,pXiFW=function(x)1,maxThreads=0,verbose=TRUE){
 DEG <- BWS <- FWS <- TRUE
	if(is.null(kit)) DEG <- FALSE
	if(!is.null(xi) && xi==0) BWS = FALSE
	if(!is.null(xiFW) && xiFW==0) {
	  FWS = FALSE
	} else if(!BWS) {
	  if(verbose) print("Backward stutter was switched on again..")
	  BWS = TRUE #Back stutter must be ON if switched off (limitation of EFMfast)
	}
	mlefit = calcMLE(nC,samples,popFreq, refData, condOrder, knownRef, kit, DEG,BWS,FWS,threshT,prC,lambda,fst,knownRel,ibd,
	                 priorBWS=pXi,priorFWS=pXiFW, maxThreads=maxThreads, verbose=verbose)
	#lims =  getParamLimits(mlefit) #obtain values under Hd if not provided otherwise
	int = calcINT(mlefit, lower=lower, upper=upper, reltol=reltol,scale=scale,maxEval=maxEval,verbose=verbose)
	return(int)
}

