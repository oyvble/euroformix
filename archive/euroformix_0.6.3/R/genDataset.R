#' @title genDataset
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description Function for generating replicated mixture samples given same contributors and model parameters.
#' @details genDataset samples random mixture peak heights given as gamma(rho*sum(h_k),tau), with h_k as peak height of k-te contributor.
#' genData conditions on alleles given by refData. Empty references are generated with population frequencies.
#' @param nC Number of contributors in model.
#' @param popFreq A list of allele frequencies for a given population.
#' @param mu Expected peak heights for a het. single contributor allele
#' @param sigma Coeffecient of variance of peak heights.
#' @param sorted Boolean for wheter sorting the contributors with respect to decreasingly mixture proportions.
#' @param threshT Required allele peak height in mixture
#' @param refData A list of given reference profiles given as refData[[i]][[s]]. Default is random from population. 
#' @param mx A vector of known mixture proportions. Default is random uniform.
#' @param nrep Number of peak height replicates (same contributors) to generate. Default is 1.
#' @param stutt A numerical stutter ratio (one ahead). Default is 0.
#' @param prC A numerical dropin probability. Default is 0.
#' @param lambda The rate parameter in the exponential distribution for simulating drop-in peak heights. Default is 0.
#' @param beta The degradation slope parameter used for simulating degradation trend (requires valid kit to be specified). Default is 1.
#' @param kit shortname of kit: {"ESX17","ESI17","ESI17Fast","ESX17Fast","Y23","Identifiler","NGM","ESSPlex","ESSplexSE","NGMSElect","SGMPlus","ESX16", "Fusion","GlobalFiler"}
#' @return List with elements theta,samples,refData where theta is the true parameters of the model. samples is a list with samples which for each samples has locus-list elements with list elements adata and hdata
#' @export

genDataset = function(nC,popFreq,mu=1000,sigma=0.1,sorted=FALSE,threshT=50,refData=NULL,mx=NULL,nrep=1,stutt=0,prC=0,lambda=0,beta=1,kit=NULL) {
  if(threshT < 0) stop("Threshold must be positive!")
  if(mu<0 || beta<0|| sigma<0) stop("Parameters must be positive!")
  if(nC<1 || round(nC)!=nC) stop("Number of contributors must be a positive whole number!")
  if(stutt<0 || stutt>1) stop("Stutter ratio must be between 0 and 1!")
  if(prC<0 || prC>1) stop("Dropin probability must be between 0 and 1!")
  nL<-length(popFreq)
  locs <- names(popFreq)
  if(is.null(locs)) locs = paste0("locus",1:nL)
  if(is.null(mx)) {
   mx <- rgamma(nC,1)
   mx=mx/sum(mx) #rdirichlet(1,rep(1,nC))  #simulate mx for contributors
  } else {
   if( sum(mx)!=1 || any(mx < 0 | mx > 1) ) stop("mx is not a simplex")
   if(length(mx)!=nC) stop("Length of mx not equal nC!")
  }
  if(sorted)  mx  <- sort(mx,decreasing=TRUE)
  mx <- c(mx)
  if(is.null(refData )) refData <- list() 

  #get kit-info
  kitinfo <- NULL
  if(!is.null(kit)) kitinfo <- getKit(kit) #get kitinfo
  if(length(kitinfo)==1) {
   print("Wrong kit specified. No degradation will be applied.")
  }

  #convert to gamma-parameters
  rho <- 1/(sigma^2)
  tau <- mu/rho

  #for each replicates:
  samples <- list()
  for(r in 1:nrep) {
   mixData <- list()
   nDropout <- nDropin <- nStutter <- rep(NA,nL) #counts dropout/dropin/stutters

   for(loc in locs) { 
#loc=locs[16]
    if( is.null(refData[[loc]])) refData[[loc]] <- list()     
    nR <- length(refData[[loc]])
    mixA <- unlist(refData[[loc]]) #vectorize
    if(nR<nC) { #sample more alleles if missing
     ran <- (nR+1):nC #new range
     for(s in ran) {
      Asim <- refData[[loc]][[s]] <-  sample(names(popFreq[[loc]]),size=2,prob=popFreq[[loc]],replace=TRUE)
      mixA = c(mixA,Asim) #keep not droped
     }
     names(refData[[loc]])[ran] <- paste0("genref",ran)
    }
    mixH = c(t(replicate(2,mx)))
    agg=aggregate(mixH,by=list(mixA),sum) #aggregate contributions
    mixA <- agg[,1]
    nA <- length(mixA) #number of alleles
    degscale <- 1 #default is no degradation scale
    if(length(kitinfo)>1) {
     locind <- toupper(kitinfo$Marker)==loc
     bp <- numeric()
     for(allel in mixA) { #for each allele
      bind <- which(locind & kitinfo$Allele==allel )
      if(length(bind)==0) { #if missing alleles in kitinfo
       close <- which.min(abs(as.numeric(allel) - as.numeric(kitinfo$Allele[locind])))
       bind <- which(locind)[close]
      } 
      bp <- c(bp,kitinfo$Size[bind]) #get fragment length
     }
     degscale = beta^((bp-125)/100) #note: must be same as in prepareC-function
    }
    mixH <- rgamma(length(agg$x),shape=rho*agg$x*degscale,scale=tau) #shape/scale given. Must be an integer!
    if(stutt>0) { #include stutter
     mixA2 <- c(agg[,1],(as.numeric(agg[,1])-1)) #stutter positions
     mixH2 <- c((1-stutt)*mixH,stutt*mixH)
     agg2=aggregate(mixH2,by=list(mixA2),sum) #aggregate stutter peaks
     mixA <- agg2[,1]
     mixH <- agg2$x
    }
    pos <- ceiling(log(1e-16)/log(prC))
    if(length(pos)>1) { #dropin done last
     prCvec <- c(1-prC/(1-prC),prC^(1:pos))
     nDI <- sample(0:pos,size=1,prob=prCvec,replace=TRUE) #number of dropin
     if(nDI>0) {
      DIA <- sample(names(popFreq[[loc]]),size=nDI,prob=popFreq[[loc]],replace=TRUE) #random from population
      if(lambda<=0) stop("Lambda was not specified greater than zero")
      DIH <- rexp(nDI,rate=lambda) + threshT #peak height of same alleles are accumulated
      mixH2 <- c(mixH,DIH)
      mixA2 <- c(mixA,DIA)
      agg2=aggregate(mixH2,by=list(mixA2),sum) #aggregate dropped in alleles
      mixA <- agg2[,1]
      mixH <- agg2$x
     }
    }
    dropped <- mixH<threshT 
    mixData[[loc]] <- list(adata=mixA[!dropped],hdata=round(mixH[!dropped]))
   } #end each locus
   names(mixData) <- locs
   samples[[r]] <- mixData
  } #end for each rep
  names(samples) <- paste0("sample",1:nrep)
  names(refData) <- locs
  return(list(theta=list(mx=mx,rho=rho,tau=tau,beta=beta),samples=samples,refData=refData))
}


