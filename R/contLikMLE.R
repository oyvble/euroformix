
#' @title contLikMLE
#' @author Oyvind Bleka
#' @description contLikMLE optimizes the likelihood function of the DNA mixture model 
#' @details Function calls the likelihood function implemented in C++ which uses the package Boost and paralellisation using OpenMP
#'
#' @param nC Number of contributors in model. Must be a constant.
#' @param samples A List with samples which for each samples has locus-list elements with list elements adata and hdata. 'adata' is a qualitative (allele) data vector and 'hdata' is a quantitative (peak heights) data vector.
#' @param popFreq A list of allele frequencies for a given population.
#' @param refData Reference objects with (2-size) allele-vector given in list element [[i]][[s]].
#' @param condOrder Specify conditioning references from refData (must be consistent order). For instance condOrder=(0,2,1,0) means that we restrict the model such that Ref2 and Ref3 are respectively conditioned as 2. contributor and 1. contributor in the model. 
#' @param knownRef Specify known non-contributing references from refData (index). For instance knownRef=(1,2) means that reference 1 and 2 is known non-contributor in the hypothesis. This affectes coancestry correction.
#' @param xi A numeric giving stutter-ratio if it is known. Default is 0, meaning stutter is not used.
#' @param prC A numeric for allele drop-in probability. Can be a vector (must contain the marker names). Default is 0.
#' @param nDone Number of optimizations required providing equivalent results (same logLik value obtained)
#' @param threshT The detection threshold given. Used when considering probability of allele drop-outs. Can be a vector (must contain the marker names).
#' @param fst The co-ancestry coefficient. Can be a vector (must contain the marker names). Default is 0.
#' @param lambda Parameter in modeled peak height shifted exponential model. Can be a vector (must contain the marker names). Default is 0.
#' @param pXi Prior function for xi-parameter (stutter). Flat prior on [0,1] is default.
#' @param delta Scaling of variation of normal distribution when drawing random startpoints. Default is 1.
#' @param kit Used to model degradation. Must be one of the shortnames of kit: check getKit()
#' @param verbose Boolean whether printing optimization progress. Default is TRUE.
#' @param maxIter Maximum number of iterations for the optimization. Default is 100.
#' @param knownRel gives the index of the reference which the 1st unknown is related to.
#' @param ibd the identical by decent coefficients of the relationship (specifies the type of relationship)
#' @param xiFW A numeric giving FW stutter-ratio if it is known.Default is 0, meaning stutter is not used.
#' @param pXiFW Prior function for xiFW-parameter (FW stutter). Flat prior on [0,1] is default.
#' @param seed The user can set seed if wanted
#' @param maxThreads Maximum number of threads to be executed by the parallelization
#' @param steptol Argument used in the nlm function for faster return from the optimization (tradeoff is lower accuracy).
#' @return ret A list(fit,model,nDone,delta,seed,prepareC) where fit is Maximixed likelihood elements for given model.
#' @export
#' @examples
#' @examples
#' \dontrun{
#' kit = "ESX17"
#' sep0 = .Platform$file.sep
#' AT0 = 50
#' popfn = paste(path.package("euroformix"),"FreqDatabases",paste0(kit,"_Norway.csv"),sep=sep0)
#' evidfn = paste(path.package("euroformix"),"examples",paste0(kit,"_3p.csv"),sep=sep0)
#' reffn = paste(path.package("euroformix"),"examples",paste0(kit,"_refs.csv"),sep=sep0)
#' popFreq = freqImport(popfn)[[1]] #obtain list with population frequencies
#' samples = sample_tableToList(tableReader(evidfn))
#' refData = sample_tableToList(tableReader(reffn))
#' dat = prepareData(samples,refData=refData,popFreq=popFreq,threshT=AT) #obtain data to use for analysis
#' plotEPG2(dat$samples,dat$refData,kit=kit,AT=AT0)
#' mlefit = contLikMLE(3,dat$samples,dat$popFreq,dat$refData,1:3,kit=kit,xi=NULL,prC=0.05,lambda=0.01,seed=1)
#' }


contLikMLE = function(nC,samples,popFreq,refData=NULL,condOrder=NULL,knownRef=NULL,xi=0,prC=0,nDone=2,threshT=50,fst=0,lambda=0,pXi=function(x)1,delta=1,kit=NULL,verbose=TRUE,maxIter=100,knownRel=NULL,ibd=c(1,0,0),xiFW=0,pXiFW=function(x)1,seed=NULL,maxThreads=32,steptol=1e-3){
  if(!is.null(seed)) set.seed(seed) #set seed if provided
    
  c <- prepareC(nC,samples,popFreq,refData,condOrder,knownRef,kit,knownRel,ibd,fst,incS=is.null(xi) || xi>0,incFS=is.null(xiFW) || xiFW>0)
  usedeg  <- !is.null(kit) #boolean whether modeling degradation TRUE=YES, FALSE=NO
  
  #Prepare fixed params:
  nM = c$nM #number of markers to evaluate

  #Check AT:
  if(length(threshT)==1) {
    ATv = rep(threshT,nM) #common detection threshold for all markers
  } else {
    ATv = setVecRightOrder(threshT, c$locs) #get right order of vector 
  }
        
  #Check dropin prob
  if(length(prC)==1) {
    pCv = rep(prC,nM) #common dropin prob for all markers
  } else {
    pCv = setVecRightOrder(prC,  c$locs) 
  }
  
  #Check hyperparam DI
  if(length(lambda)==1) {
    lambdav = rep(lambda,nM) #common hyperparam DI for all markers
  } else {
    lambdav = setVecRightOrder(lambda,  c$locs) 
  }
  
  #Check theta correction
  if(length(fst)==1) {
    fstv = rep(fst,nM) #common theta correction for all markers
  } else {
    fstv = setVecRightOrder(fst,  c$locs)
  }
  names(ATv) <- names(pCv) <- names(lambdav) <- names(fstv) <- c$locs #insert locus names (correct order)
  
  #Prefitting data based on the model for sum of the peak heights  to find proper startvalues for MLE 
  sumY <- meanbp <- rep(NA,nM)
  startind1 = 0 #start index (no reps)
  startindR = 0 #start index (all reps)
  for(i in 1:nM) { #for each marker
    nAreps = c$nA[i]*c$nReps[i] #number of alleles over all reps
    ind1 = startind1 + c(1:c$nA[i]) #obtain index of PH vector to use
    indR = startindR + c(1:nAreps) #obtain index of PH vector to use
    sumY[i] <- sum(c$YvecLong[indR])/c$nReps[i] #take sum of the peak heights (average over number of reps)
    meanbp[i] <- mean(c$basepairLong[ind1]) #take sum of the peak heights
    startind1 =  tail(ind1,1)
    startindR =  tail(indR,1)
  }
  if(usedeg) { #if degradation 
    alpha0 <- mean(sapply(samples,function(x) max(sapply(x, function(y) sum(as.numeric(y$hdata)))))/(2)) #mean het. allele. peak height at largest marker when degradation
    th0 <- c(alpha0,0.4,0.8)
    suppressWarnings({ tryCatch({  th0 <- euroformix::fitgammamodel(sumY,x=meanbp,niter=1,delta=0,offset=0,scale=1)  }, error = function(e) e) }) #fit model with DEG
  } else { #if no degradation
    alpha0 <- mean(sapply(samples,function(x) mean(sapply(x, function(y) sum(as.numeric(y$hdata)))))/(2)) #mean het. allele. peak height averaged on all markers used when no degradation
    th0 <- c(alpha0,0.4)
    suppressWarnings({ tryCatch({ th0 <- euroformix::fitgammamodel(sumY,niter=1,delta=0)  }, error = function(e) e) }) #fit sum of peak height model
  }
  #if(verbose) cat(paste0("theta0=",paste0(th0,collapse=",")))  
  
  #PREPEARING THE LIKELIHOOD OPTIMZATION (what parameters are provided?):
  useParamOther = rep(TRUE,3)  #index of parameters used (set NA if fixed)
  if(!usedeg) useParamOther[1] = FALSE #degrad not used
  if(!is.null(xi)) useParamOther[2] = FALSE #BW stutter not used
  if(!is.null(xiFW)) useParamOther[3] = FALSE  #FW stutter not used
  indexParamOther = rep(NA,3)
  if(any(useParamOther)) indexParamOther[useParamOther] = 1:sum(useParamOther) #init indices
  
  #function for calling on C-function: Must convert "real domain" values back to model params 
  negloglikYphi <- function(phi,progressbar=TRUE) { #assumed order: mixprop(1:C-1),mu,sigma,beta,xi
    mixprop = 1 #default variable (not used in c code if nC=1)
    if(nC>1) {
      mixprop = phi[1:(nC-1)]#mixture proportions are converted back in the c code
    } #end if more than 1 contr
    mu1 = exp(phi[nC])
    sig1 = exp(phi[nC+1])
    
    beta1 = 1.0 #set default values
    xiB = xi #set default values
    xiF = xiFW #set default values
    tmp =  1/(1+exp(-phi[(nC+1) + indexParamOther ])) #extract other params 
    if(useParamOther[1]) beta1 = tmp[1]
    if(useParamOther[2]) xiB = tmp[2]
    if(useParamOther[3]) xiF = tmp[3]
    
    muv  = rep(mu1,nM)
    sigmav =  as.numeric(rep(sig1,nM))
    betav =  rep(beta1,nM) 
    xiBv = rep(xiB,nM) 
    xiFv = rep(xiF,nM) 
    
    loglik = .C("loglikgammaC",as.numeric(0),c$nC,c$NOK,c$knownGind,as.numeric(mixprop),as.numeric(muv),as.numeric(sigmav),as.numeric(betav),as.numeric(xiBv),as.numeric(xiFv),as.numeric(ATv),as.numeric(pCv),as.numeric(lambdav),as.numeric(fstv),c$nReps,c$nM,c$nA,c$YvecLong,c$FvecLong,c$nTypedLong,c$maTypedLong,c$basepairLong,c$BWvecLong,c$FWvecLong,c$nPS,c$BWPvecLong,c$FWPvecLong,as.integer(maxThreads),as.integer(1),c$anyRel,c$relGind,c$ibdLong,PACKAGE="euroformix")[[1]]
    if(is.null(xi))  loglik <- loglik + log(pXi(xiB)) #weight with prior of xi
    if(is.null(xiFW))  loglik <- loglik + log(pXiFW(xiF)) #weight with prior of xiFW
    
    if(verbose && progressbar) { #only show progressbar if verbose(and if decided to show)
      progcount <<- progcount + 1 #update counter
      #tcltk::setTkProgressBar(progbar,progcount,label=paste0(expUpperTime,"\n", round(progcount/maxIterProgress*100, 0),"% done"))
      setTxtProgressBar(progbar,progcount)
    } 
    return(-loglik) #weight with prior of stutter.
  }
  
  np1 = nC - 1 #number of mixture params
  np2 = sum(usedeg) + sum(is.null(xi)) + sum(is.null(xiFW)) #number of beta/xi params
  np <- np1 + np2 + 2 #number of unknown parameters
  
  ncond = sum(condOrder>0) #number of conditional
  nU = nC - ncond #number of unknows
  
  logit = function(x) log(x/(1-x))
  paramrandomizer = function(verbose2=FALSE) { #draw "good" randoms start points (good guess)
    
    #CONSIDER Mixture porpotions    
    mxrnd = rgamma(nC,1) #Draw simplex (flat)
    mxrnd = mxrnd/sum(mxrnd)
    if(nU>1) { #sort if more than 1 unknown
      ind = (ncond+1):nC #sort Mix-prop for the unknowns 
      mxrnd[ind] = sort(mxrnd[ind],decreasing = TRUE) #sort Mx in decreasing order
    }
    #mxrnd = mxrnd[-nC] #don't use mixprop for last ind (not needed!)
    randParam = mxrnd #obtain random params (used for printing)
    
    #convert Mx values to real domain (nu:
    nurnd = numeric()
    if(nC>1) {
      cs = 0 #c( 0,cumsum(mxrnd)) #cumulative sum of mixture proportins
      for(cc in 1:(nC-1)) { #traverse contributors (Restricted)
        nurnd = c(nurnd, logit( mxrnd[cc]/(1-cs))) 
        cs = cs + mxrnd[cc] #update sum
      }
    }
    
    #CONSIDER PH prop variables
    th1 = th0[1:2]  #PH prop variables
    sdPH = delta*0.15*th1 #obtain considered SD of PH props 
    PHrnd = abs( rnorm(2,th1,sd=sdPH))
    logPHrnd = log(PHrnd)  #Obtain random start for mu/sigma, Note using the delta here (should be small)
    randParam = c(randParam,PHrnd) #add random for PHprop
    
    #CONSIDER other variables (degrad/BW/FW)
    otherrnd = numeric() #random for beta,xi etc
    if(usedeg) {
      maxval = 5 #maximum transformed degrad slope param
      degval = logit( th0[3] ) #get transformed degrad slope param value
      if( is.infinite(degval) || degval>maxval ) degval = maxval #insert fixed large value (maxval=5)
      otherrnd = c(otherrnd, degval ) #extract for degradation
    }
    if(is.null(xi)) otherrnd = c(otherrnd, logit(0.05) ) #assume stutter prop expectation of 0.05
    if(is.null(xiFW)) otherrnd = c(otherrnd, logit(0.01) ) #assume stutter prop expectation of 0.025
    
    if(length(otherrnd)>0) otherrnd = rnorm(length(otherrnd),otherrnd,sd=0.5) #Note small sd (because shifted with expected trend)
    randParam = c(randParam, 1/(1+exp(-otherrnd))) #transform back and add
    
    if(verbose2)   cat(paste0("\nRandom startparam=",paste0(prettyNum(randParam),collapse="/")))
    
    return(c(nurnd,logPHrnd,otherrnd))    
  }
  
  secondToTimeformat <- function(t){ #converts seconds to time format
    paste(formatC(t %/% (60*60) %% 24, width = 2, format = "d", flag = "0"), #hours
          formatC(t %/% 60 %% 60, width = 2, format = "d", flag = "0"), #mins
          formatC(t %% 60, width = 2, format = "d", flag = "0"),sep = ":") #seconds
   }
  
  #Strategy to obtain global maximum: Assure that maximum logLik (maxL) is obtained in 'nDone' optimizations 
  maxITERS <- 30 #number of possible times to be INF or not valid optimum before any acceptance
  maxL <- -Inf #value of maximum obtained loglik
  nOK <- 0 #number of times for reaching largest previously seen optimum
  nITER <- 0 #number of times beeing INF (invalid value)
  
  iterscale = np-2 #get iteration scale in MLE (number of params - 2) 
  maxIterProgress = iterscale*maxIter #maximum iterations in progressbar
  
  suppressWarnings({
    while(nOK<nDone) {
      #Obtain random start values for parameters
      
      p0 <- paramrandomizer(verbose2=FALSE) #generate random start value on Real (don't need to)
      
      timeOneCall = system.time({ #estimate the time for calling the likelihood function one time 
        likval <- negloglikYphi(phi=p0,FALSE)   #check if start value was accepted
      })[3] #obtain time in seconds

      if( is.infinite(likval) ) { #if it was infinite (invalid)
        nITER = nITER + 1	 
      } else {
        if(verbose) {
          showProgressBar = FALSE #show progress bar only if calculatations take more than 10 seconds
          expectedTimeProgress0 = timeOneCall*maxIterProgress
          if(expectedTimeProgress0 > 10) showProgressBar = TRUE #show progress if upper time >10s
          if(showProgressBar) {
            expectedTimeProgress = secondToTimeformat(expectedTimeProgress0)
            cat(paste0("\nExpected (upper) time is ", expectedTimeProgress, " (HH:MM:SS):\n")) 
            #\n---------------------------------------------------
          }
        }
        
        progcount  = 1 #progress counter
        #if(verbose) progbar <- tcltk::tkProgressBar(min = 0, max = iterscale*maxIter,width = 300) #create progress bar
        if(verbose && showProgressBar) progbar <- txtProgressBar(min = 0, max = maxIterProgress, style = 3) #create progress bar
        
        tryCatch( {
          foo <- nlm(f=negloglikYphi, p=p0,hessian=TRUE,iterlim=maxIter,steptol=steptol, progressbar=showProgressBar)#,print.level=2)
          Sigma <- solve(foo$hessian)
          
          if(all(diag(Sigma)>=0) && foo$iterations>2) { #} && foo$code%in%c(1,2)) { #REQUIREMENT FOR BEING ACCEPTED
            nITER <- 0 #reset INF if accepted
            likval <- -foo$min #obtain local maximum
            
            #was the maximum (approx) equal the prev: Using decimal numbers as difference 
            isEqual = !is.infinite(maxL) && abs(likval-maxL)<0.01 #all.equal(likval,maxL, tolerance = 1e-2) # # #was the maximum (approx) equal the prev?
            
            if(isEqual) { 
              if(verbose) {
                if(showProgressBar) print("") #skip line
                print(paste0("Equal maximum found: loglik=",likval))
              }
              nOK = nOK + 1 #add counter by 1
            } else {  #if values were different
              if(likval>maxL) { #if new value is better
                nOK = 1 #first accepted optimization found
                maxL <- likval #maximized likelihood
                maxPhi <- foo$est #set as topfoo     
                maxSigma <- Sigma 
                if(verbose) {
                  if(showProgressBar) print("") #skip line
                  print(paste0("New maximum at loglik=",likval))
                } 
                
                # if(verbose) cat(paste0("maxPhi=",paste0(maxPhi,collapse=","))) #removed in v2.0
              } else {
                if(verbose) {
                  if(showProgressBar) print("") #skip line
                  print(paste0("Local (non-global) maximum found at logLik=",likval))
                }
              }
            } 
            if(verbose) print(paste0(" (",nOK,"/",nDone,") optimizations done"))
            #flush.console()
          } else { #NOT ACCEPTED
            nITER <- nITER + 1 
          }
        },error=function(e) e,finally = {nITER <- nITER + 1} ) #end trycatch (update counter)
        #if(verbose) close(progbar) #NECESSARY FOR TCLTK
        
      } #end if else 
      
      if(nOK==0 && nITER>maxITERS) {
        nOK <- nDone #finish loop
        maxL <- -Inf #maximized likelihood
        maxPhi <- rep(NA,np) #Set as NA
        maxSigma <- matrix(NA,np,np)#Set as NA
        break #stop loop if too many iterations  
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
  otherInd = setdiff(1:np,(1:(nC+1))) #get indices of other params
  
  if( length(otherInd)>0  ) { #if remaining params: inverse-logit of remaining params
    tmp = 1/(1+exp(-maxPhi[otherInd]))
    thetahat = c(thetahat , tmp)
    thetahat2 = c(thetahat2 , tmp)
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
    
    if(  length(otherInd)>0   ) { #if remaining params: 
      tmp = exp(-phi[otherInd])
      J[cbind(otherInd,otherInd)] = tmp*(1+tmp)^(-2) #insert values (vectorized)
    }
    
    return(J)
  } #end jacobian
  J <- Jacob(phi=maxPhi,mle=thetahat)
  Sigma <- (t(J)%*%maxSigma%*%J) #this is correct covariance of thetahat. Observed hessian is used
  
  #get extended Sigma (all parameters)
  Sigma2 <- matrix(NA,nrow=np+1,ncol=np+1) #extended covariance matrix also including mx[nC]
  Sigma2[nC:np+1,nC:np+1] <- Sigma[nC:np,nC:np] 
  Sigma2[nC:np+1,nC:np+1] <- Sigma[nC:np,nC:np] 
  if(nC>1) {
    Sigma2[nC:np+1,1:(nC-1)] <- Sigma[nC:np,1:(nC-1)] 
    Sigma2[1:(nC-1),1:(nC-1)] <- Sigma[1:(nC-1),1:(nC-1)] 
    Sigma2[1:(nC-1),nC:np+1] <- Sigma[1:(nC-1),nC:np] 
    Sigma2[nC,nC] <- sum(Sigma[1:(nC-1),1:(nC-1)])
    for(k in (1:(np+1))[-nC]) {
      Sigma2[nC,k] <- Sigma2[k,nC] <- -sum(Sigma[1:(nC-1),k-sum(k>nC)]) 
    }
  } else {
    Sigma2[1,1:(np+1)] <- Sigma2[1:(np+1),1] <- 0 #no uncertainty
  }
  
  #Standard error for theta:
  thetaSE <- sqrt(diag(Sigma2))
  
  mxName <- "Mix-prop. C" #"mx"
  muName <- "P.H.expectation" #"mu"
  sigmaName <- "P.H.variability" #"sigma"
  betaName <- "Degrad. slope"#"beta"
  xiName <- "BWstutt-prop."#"xi"
  xiFWName <- "FWstutt-prop."#"xi"
  
  thetanames0 <- c(muName,sigmaName)
  phinames <- paste0("log(",thetanames0,")")
  if(nC>1) {
    phinames  <- c(paste0("nu",1:(nC-1)),phinames)
    thetanames <- c(paste0(mxName,1:(nC-1)),thetanames0)
  } else {
    thetanames <- thetanames0
  }
  thetanames2 <- c(paste0(mxName,1:nC),thetanames0)
  
  if( length(otherInd)>0 ) { #if other indices
    othernames = character()
    if(useParamOther[1]) othernames = c(othernames,betaName)
    if(useParamOther[2]) othernames = c(othernames,xiName)
    if(useParamOther[3]) othernames = c(othernames,xiFWName)
    
    phinames <- c(phinames,paste0("logit(",othernames,")") )
    thetanames <- c(thetanames, othernames )
    thetanames2 <- c(thetanames2, othernames )
    
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
  nU <- nC #-ret$nK #number of unknowns
  if(nU>1) { #if more than 1 unknown 
    logmargL <- log(factorial(nU)) + logmargL #get correct ML adjusting for symmetry
  }
  fit <- list(phihat=maxPhi,thetahat=thetahat,thetahat2=thetahat2,phiSigma=maxSigma,thetaSigma=Sigma,thetaSigma2=Sigma2,loglik=maxL,thetaSE=thetaSE,logmargL=logmargL)
  #store model:
  model <- list(nC=nC,samples=samples,popFreq=popFreq,refData=refData,condOrder=condOrder,knownRef=knownRef,xi=xi,prC=pCv,threshT=ATv,fst=fstv,lambda=lambdav,pXi=pXi,kit=kit,knownRel=knownRel,ibd=ibd,pXiFW=pXiFW,xiFW=xiFW)
  ret <- list(fit=fit,model=model,nDone=nDone,delta=delta,steptol=steptol,seed=seed,prepareC=c) #store seed
  return(ret)
} #end function

