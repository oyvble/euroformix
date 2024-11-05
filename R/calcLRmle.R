#' @title calcLRmle
#' @description A function for summarizing LR (MLE)
#' @param mlefitHp A fitted object returned from calcMLE (under Hp)
#' @param mlefitHd A fitted object returned from calcMLE (under Hd)
#' @export

calcLRmle = function(mlefitHp, mlefitHd) {
  lr = (mlefitHp$fit$loglik-mlefitHd$fit$loglik)
  lri = c(mlefitHp$fit$logliki-mlefitHd$fit$logliki)
  markersHp = names(mlefitHp$fit$logliki)
  markersHd = names(mlefitHd$fit$logliki)
  if(!all(markersHp==markersHd)) warning("Note: Marker names does not coincide for the LR-per marker!")
  lri = setNames(lri,markersHd)
  
  #CALCULATE OTHER LR based statistics (alternatives)
  LRlap <- exp(mlefitHp$fit$logmargL - mlefitHd$fit$logmargL)#/log(10) #calculate laplace approximated LRs
  
  #Calculated Adjusted LR (sub source) based on number of unknowns with unequal Mx:
  MxHp = mlefitHp$fit$thetahat2[1: mlefitHp$model$nC ] #mixture prop est Hp
  MxHd = mlefitHd$fit$thetahat2[1: mlefitHd$model$nC ] #mixture prop est Hd
  MxUp = tail(MxHp,mlefitHp$model$nU) #get Mx for Unknowns Hp
  MxUd = tail(MxHd,mlefitHd$model$nU) #get Mx for Unknowns Hd
  epsround = 4 #number of decimals used to be equal in Mx
  nUpDiffMx = length(unique(round(MxUp,epsround))) #get length of unique Mx of unknowns (hp)
  nUdDiffMx = length(unique(round(MxUd,epsround))) #get length of unique Mx of unknowns (hd)
  LRadj = exp(lr + lgamma(nUpDiffMx+1) - lgamma(nUdDiffMx+1) ) #obtain adjusted lr
  
  #Calculate upper theoretical LR limit for POI
  condHp = mlefitHp$model$condOrder 
  condHpIdx = which(mlefitHp$model$condOrder>0)
  condHdIdx = which(mlefitHd$model$condOrder>0)
  POIidx = condHp[setdiff(condHpIdx,condHdIdx)] #get ref-ind in hp but not in hd
  LRupper = euroformix::calcLRupper(POIidx,mlefitHp) #obtain calculated LR. NOTE: NEED TO SEND hp fit AND position of POIidx
  
  ret = list(LR=exp(lr),LRmarker=exp(lri),log10LR=lr/log(10),log10LRmarker=lri/log(10),LRlap=LRlap, LRadj=LRadj,
	LRupper=10^LRupper, log10LRupper=LRupper)
  return(ret)
}