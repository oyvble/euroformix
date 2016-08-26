#' @title logLiki
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description logLiki calculates the likelihood of each marker of the STR DNA mixture given a theta
#' @details
#' Function calls procedure in c++ by using the package Armadillo and Boost.
#'
#' @param mlefit Fitted object using contLikMLE
#' @return ret A vector with log-likelihood-values for each locus for given model
#' @export

#Get value of likelihood for each marker given theta
logLiki <- function(mlefit){
 theta <- mlefit$fit$thetahat #condition on mle parameter
 model <- mlefit$model #take out assumed model with given data
 np <- length(theta) #number of unknown parameters
 nC <- model$nC
 locs <- names(model$popFreq)
 nL <- length(locs)
 xi <- model$xi
 nodeg <- is.null(model$kit) #check for degradation

 theta2 <- theta[1:(nC+1)] #take out mx,mu,sigma
 if(nodeg) {
   theta2 <- c(theta2,1) #add beta=1 to parameters. Means no degradation
 } else {
   theta2 <- c(theta2,theta[nC+2]) #add beta to parameters
 }
 if(is.null(xi)) {  #if xi unknown
   theta2 <- c(theta2,theta[np]) #add xi param to parameters
 } else { #if xi known
   theta2 <- c(theta2,as.numeric(xi)) #add xi param to parameters
 }

 if(any(is.na(theta))) {
  logLi <- rep(NA,nL)
 } else {
  logLi <- numeric()
  logliktheta <- function() {   #call c++- function: 
    Cval  <- .C("loglikgammaC",as.numeric(0),as.numeric(theta2),as.integer(np),ret$nC,ret$nK,ret$nL,ret$nS,ret$nA,ret$obsY,ret$obsA,ret$CnA,ret$allAbpind,ret$nAall,ret$CnAall,ret$Gvec,ret$nG,ret$CnG,ret$CnG2,ret$pG,ret$pA, as.numeric(model$prC), ret$condRef,as.numeric(model$threshT),as.numeric(model$fst),ret$mkvec,ret$nkval,as.numeric(model$lambda),ret$bp,as.integer(0),PACKAGE="euroformix")[[1]]
    if(is.null(xi)) Cval <- Cval + log(model$pXi(theta[np]))
    return(Cval) 
  }
  for(loc in locs) { #traverse for each locus
   samples <- lapply(model$samples,function(x) x[loc])
   ret <- prepareC(nC=nC,samples,popFreq=model$popFreq[loc],refData=model$refData[loc],condOrder=model$condOrder,knownRef=model$knownRef,kit=model$kit)
   logLi <- c(logLi,logliktheta())
  }
 } #end case
 names(logLi) <- locs
 return(logLi)
}



