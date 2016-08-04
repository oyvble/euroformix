
#' @title contLikMLE
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description contLikMLE optimizes the likelihood of the STR DNA mixture given some assumed a bayesian model.
#' @details The procedure are doing numerical optimization to approximate the marginal probabilit over noisance parameters. Mixture proportions have flat prior.
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
#' @return ret A list(fit,model,nDone,delta) where fit is Maximixed likelihood elements for given model.
#' @export
#' @references Cowell,R.G. et.al. (2014). Analysis of forensic DNA mixtures with artefacts. Applied Statistics, 64(1),1-32.
#' @keywords continuous model, Maximum Likelihood Estimation
contLikMLE = function(nC,samples,popFreq,refData=NULL,condOrder=NULL,knownRef=NULL,xi=NULL,prC=0,nDone=1,threshT=50,fst=0,lambda=0,pXi=function(x)1,delta=10,kit=NULL,verbose=TRUE){
 ret <- prepareC(nC,samples,popFreq,refData,condOrder,knownRef,kit)
 nodeg  <- is.null(kit) #boolean whether modeling degradation FALSE=YES, TRUE=NO

 #Prefitting data based on the model for sum of the peak heights  to find proper startvalues for MLE 
 sumY <- meanbp <- rep(NA,ret$nL*ret$nS)
 for(i in 1:ret$nL) {
  for(j in 1:ret$nS) {
   ind <- (ret$nS*(i-1) + j)
   rng1 <- (ret$CnA[ind]+1):ret$CnA[ind+1]
   rng2 <- (ret$CnAall[i]+1):ret$CnAall[i+1]
   sumY[ind] <- sum(ret$obsY[rng1]) #take sum of the peak heights
   meanbp[ind] <- mean(ret$bp[rng2]) #take sum of the peak heights
  }
 }
# plot(meanbp,sumY)
 negloglik <- function(th) {
  th <- exp(th)
  if(!nodeg)  val <- -sum(dgamma(sumY,shape=2/th[2]^2*th[3]^(meanbp),scale=th[1]*th[2]^2,log=TRUE)) 
  if(nodeg)  val <- -sum(dgamma(sumY,shape=2/th[2]^2,scale=th[1]*th[2]^2,log=TRUE)) 
  if(is.infinite(val)) val <- .Machine$integer.max 
  return(val)
 }
 foo <- NULL
 if(!nodeg) suppressWarnings({ tryCatch({ foo <- nlm(negloglik, c(log(mean(sumY)),0.3,0.8) ) }, error = function(e) e) })
 if(nodeg) suppressWarnings({ tryCatch({ foo <- nlm(negloglik, c(log(mean(sumY)),0.3) ) }, error = function(e) e) })
 if(!is.null(foo)) {
  th0 <- exp(foo$est)
 } else {
  if(nodeg) alpha0 <- mean(sapply(samples,function(x) mean(sapply(x, function(y) sum(as.numeric(y$hdata)))))/(2)) #mean het. allele. peak height averaged on all markers used when no degradation
  if(!nodeg) alpha0 <- mean(sapply(samples,function(x) max(sapply(x, function(y) sum(as.numeric(y$hdata)))))/(2)) #mean het. allele. peak height at largest marker when degradation
  th0 <-  c(alpha0,0.4) #guess sigma param
  if(!nodeg)   th0 <-  c(th0,0.8) #guess beta param 
 }
 if(verbose) print(paste0("theta0=",paste0(th0,collapse=",")))  

 #function for calling on C-function:
 negloglikYphi <- function(phi) {   
  phi2 <- phi[1:(nC+1)] #take out mx,mu,sigma
  if(nodeg) {
    phi2 <- c(phi2,0) #add beta=0 to parameters
  } else {
    phi2 <- c(phi2,phi[nC+2]) #add beta to parameters
  }
  if(is.null(xi)) {  #if xi uknown
    phi2 <- c(phi2,phi[np2]) #add xi param to parameters
  } else { #if xi known
    phi2 <- c(phi2,xi) #add xi param to parameters
  }
  loglik <- .C("loglikgammaC",as.numeric(0),as.numeric(phi2),as.integer(np),ret$nC,ret$nK,ret$nL,ret$nS,ret$nA,ret$obsY,ret$obsA,ret$CnA,ret$allAbpind,ret$nAall,ret$CnAall,ret$Gvec,ret$nG,ret$CnG,ret$CnG2,ret$pG,ret$pA, as.numeric(prC), ret$condRef,as.numeric(threshT),as.numeric(fst),ret$mkvec,ret$nkval,as.numeric(lambda),as.numeric(ret$bp),as.integer(1),PACKAGE="euroformix")[[1]]
  if(is.null(xi))  loglik <- loglik + log(pXi(1/(1+exp(-phi[np2])))) #weight with prior of tau and 
  return(-loglik) #weight with prior of tau and stutter.
 }
 maxITER <- 100 #number of possible times to be INF or not valid optimum before any acceptance
 np2 <- np <- nC + 2 + sum(is.null(xi)) #number of unknown parameters (including degeneration param)
 if(nodeg) np2 <- np2 - 1
 maxL <- -Inf #
 nOK <- 0 #number of times for reaching optimum
 nITER <- 0 #number of times beeing INF
 suppressWarnings({
  while(nOK<nDone) {
    p0 <- rnorm(np2,sd=delta) #generate random start value on Real
    p0[nC:(nC+length(th0)-1)] <- rnorm(length(th0),log(th0),sd=delta) #generate random start value for mu with mean alpha0
#    p0[nC] <- rnorm(1,log(alpha0),sd=delta) #generate random start value for mu with mean alpha0
    likval <- negloglikYphi(p0)  
    if(is.infinite(likval)) { #if it was infinite
	 nITER <- nITER + 1	 
    } else {
     tryCatch( {
       foo <- nlm(f=negloglikYphi, p=p0,hessian=TRUE)
       Sigma <- solve(foo$hessian)
       if(all(diag(Sigma)>=0) && foo$code%in%c(1,2) && foo$iterations>2) { #REQUIREMENT FOR BEING ACCEPTED
     	  nITER <- 0 #reset INF if accepted
        likval <- -foo$min
        nOK=nOK+1 #it was accepted as an optimum
        if(likval>maxL) {
         maxL <- likval #maximized likelihood
         maxPhi <- foo$est #set as topfoo     
         maxSigma <- Sigma 
         if(verbose) print(paste0("loglik=",maxL))
         if(verbose) print(paste0("maxPhi=",paste0(maxPhi,collapse=",")))
        }
	  if(verbose) print(paste0("Done with ",nOK,"/",nDone," optimizations"))
        flush.console()
       } else { #NOT ACCEPTED
     	  nITER <- nITER + 1 
       }
     },error=function(e) e) #end trycatch 
    }
    if(nOK==0 && nITER>maxITER) {
     nOK <- nDone #finish loop
     maxL <- -Inf #maximized likelihood
     maxPhi <- rep(NA,np2) #Set as NA
     maxSigma <- matrix(NA,np2,np2)#Set as NA
    }
  } #end while loop
 })
 
 #transfer back: 
 mx <- numeric()
 if(nC>1) {
  mx <- 1/(1+exp(-maxPhi[1:(nC-1)]))
  if(nC>2) { #need to transfer back
   for(i in 2:(nC-1)) {
    mx[i] <- mx[i]*(1-sum(mx[1:(i-1)]))
   }
  }
 }
 tmp <- exp(maxPhi[nC:(nC+1)]) #inverse-log
 thetahat <- c(mx,tmp) #last index is removed. This could again be a known contributor
 thetahat2 <- c(mx,1-sum(mx),tmp) #last index is removed. This could again be a known contributor
 if(!nodeg) {
  tmp <- exp(maxPhi[nC+2])
  thetahat <- c(thetahat,tmp) #add beta to parameters
  thetahat2 <- c(thetahat2,tmp) #add beta to parameters
 }
 if(is.null(xi)) {
  tmp <- 1/(1+exp(-maxPhi[np2])) #inverse-logit
  thetahat <- c(thetahat,tmp) #add xi to parameters
  thetahat2 <- c(thetahat2,tmp) #add xi to parameters
 }

 #Delta-method to Sigma matrix
 Jacob <- function(phi,mle) { #Jabobian matrix (in value phi)
  J <- matrix(0,length(phi),length(phi))
  if(nC>1) {
   DmDm <- matrix(0,nC-1,nC-1)
   irange <- 1:(nC-1) #range of mixture proportions
   xtmp <- 1/(1+exp(-phi[irange ])) #aux variable
   dxtmp <-  exp(-phi[irange])*xtmp^2 #derivative of aux variable
   for(i in irange) {
    for(j in 1:i) {
     if(j==i) {
       DmDm[i,i] <- dxtmp[i]
       if(i>1) DmDm[i,i] <- DmDm[i,i]*(1-sum(mle[1:(j-1)])) #note using mle(theta) here!
     } else { #cross-derivatives
       DmDm[i,j] <- -xtmp[i]*sum(DmDm[1:(i-1),j])
     }
    } #end for each col j (xj)
   } #end for each row i (fi)
  J[1:(nC-1),1:(nC-1)] <- DmDm
  }
  for(i in nC:(nC+1)) J[i,i] <- exp(phi[i])
  if(!nodeg) J[nC+2,nC+2] <- exp(phi[nC+2])
  if(is.null(xi)) {
   tmp <- exp(-phi[np2])
   J[np2,np2] <- tmp*(1+tmp)^(-2)
  }
  return(J)
 } #end jacobian
 J <- Jacob(maxPhi,mle=thetahat)
 Sigma <- (t(J)%*%maxSigma%*%J) #this is correct covariance of thetahat. Observed hessian is used
 #Sigma <- (J%*%maxSigma%*%t(J))/sqrt(np) #this formula is not correct

 #get extended Sigma (all parameters)
 Sigma2 <- matrix(NA,nrow=np2+1,ncol=np2+1) #extended covariance matrix also including mx[nC]
 Sigma2[nC:np2+1,nC:np2+1] <- Sigma[nC:np2,nC:np2] 
 Sigma2[nC:np2+1,nC:np2+1] <- Sigma[nC:np2,nC:np2] 
 if(nC>1) {
  Sigma2[nC:np2+1,1:(nC-1)] <- Sigma[nC:np2,1:(nC-1)] 
  Sigma2[1:(nC-1),1:(nC-1)] <- Sigma[1:(nC-1),1:(nC-1)] 
  Sigma2[1:(nC-1),nC:np2+1] <- Sigma[1:(nC-1),nC:np2] 
  Sigma2[nC,nC] <- sum(Sigma[1:(nC-1),1:(nC-1)])
  for(k in (1:(np2+1))[-nC]) {
   Sigma2[nC,k] <- Sigma2[k,nC] <- -sum(Sigma[1:(nC-1),k-sum(k>nC)]) 
  }
 } else {
  Sigma2[1,1:(np2+1)] <- Sigma2[1:(np2+1),1] <- 0 #no uncertainty
 }

 #Standard error for theta:
 thetaSE <- sqrt(diag(Sigma2))


 thetanames0 <- c("mu","sigma")
 if(!nodeg) thetanames0 <- c(thetanames0,"beta")
 phinames <- paste0("log(",thetanames0,")")

 if(nC>1) {
  phinames  <- c(paste0("nu",1:(nC-1)),phinames)
  thetanames <- c(paste0("mx",1:(nC-1)),thetanames0)
 } else {
  thetanames <- thetanames0
 }
 thetanames2 <- c(paste0("mx",1:nC),thetanames0)
 if(is.null(xi)) {
  phinames <- c(phinames,"logit(xi)")
  thetanames <- c(thetanames,"xi") 
  thetanames2 <- c(thetanames2,"xi") 
 }
 colnames(maxSigma) <- rownames(maxSigma) <- phinames 
 colnames(Sigma) <- rownames(Sigma) <- thetanames
 colnames(Sigma2) <- rownames(Sigma2) <- thetanames2
 names(maxPhi) <- phinames
 names(thetahat) <- thetanames
 names(thetahat2) <- thetanames2
 names(thetaSE) <- thetanames2

 #laplace approx:
 logmargL <- 0.5*(np*log(2*pi)+determinant(Sigma)$mod[1]) + maxL #get log-marginalized likelihood
 nU <- nC-ret$nK #number of unknowns
 if(nU>1) { #if more than 1 unknown 
  logmargL <- log(factorial(nU)) + logmargL #get correct ML adjusting for symmetry
 }
 fit <- list(phihat=maxPhi,thetahat=thetahat,thetahat2=thetahat2,phiSigma=maxSigma,thetaSigma=Sigma,thetaSigma2=Sigma2,loglik=maxL,thetaSE=thetaSE,logmargL=logmargL)
 #store model:
 model <- list(nC=nC,samples=samples,popFreq=popFreq,refData=refData,condOrder=condOrder,knownRef=knownRef,xi=xi,prC=prC,threshT=threshT,fst=fst,lambda=lambda,pXi=pXi,kit=kit)
 ret <- list(fit=fit,model=model,nDone=nDone,delta=delta)
 return(ret)
} #end function

