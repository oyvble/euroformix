#' @title genDataset
#' @author Oyvind Bleka 
#' @description Function for generating mixture samples for given model parameters
#' @details genDataset samples random mixture peak heights given as gamma(rho*sum(h_k),tau), with h_k as peak height of k-te contributor.
#' genData conditions on alleles given by refData. Empty references are generated with population frequencies.
#' @param nC Number of contributors in model.
#' @param popFreq A list with allele frequencies for a given population.
#' @param mu Expected peak heights for a het. single contributor allele
#' @param sigma Coeffecient of variance of peak heights.
#' @param sorted Wheter sorting the contributors with respect to decreasingly mixture proportions.
#' @param threshT Required allele peak height in mixture (can be a vector with names giving the loci names)
#' @param refData A list with given reference profiles given as refData[[i]][[s]]. Default is random from population. 
#' @param mx A vector with known mixture proportions. Default is random uniform.
#' @param nrep Number of peak height replicates (same contributors) to generate. Default is 1.
#' @param stutt A numerical stutter proportion (n-1). Default is 0.
#' @param prC A numerical dropin probability (can be a vector with names giving the loci names). Default is 0. 
#' @param lambda The rate parameter in the exponential distribution for simulating drop-in peak heights (can be a vector with names giving the loci names). Default is 0.
#' @param beta The degradation slope parameter used for simulating degradation trend (requires valid kit to be specified). Default is 1.
#' @param kit shortname of kit: "ESX17","ESI17","ESI17Fast","ESX17Fast","Y23","Identifiler","NGM","ESSPlex","ESSplexSE","NGMSElect","SGMPlus","ESX16", "Fusion","GlobalFiler"
#' @param stuttFW A numerical Forward stutter proportion (n+1). Default is 0.
#' @return List with elements theta,samples,refData where theta is the true parameters of the model. samples is a list with samples which for each samples has locus-list elements with list elements adata and hdata
#' @export
#' @examples
#' \dontrun{ 
#' kit = "SGMPlus"
#' popfile = paste(path.package("euroformix"),"FreqDatabases",
#'  paste0(kit,"_Norway.csv"),sep=.Platform$file.sep)
#' popFreq = freqImport(popfile)[[1]] #obtain list with population frequencies
#' gen = genDataset(2,popFreq,beta=.7,kit=kit,stutt=.1,stuttFW=.05)
#' plotEPG2(gen$samples,kit = kit,refData=gen$refData,AT=50) #visualize samples
#' }

genDataset = function(nC,popFreq,mu=1000,sigma=0.1,sorted=FALSE,threshT=50,refData=NULL,mx=NULL,nrep=1,stutt=0,prC=0,lambda=0,beta=1,kit=NULL,stuttFW=0) {
  if(any(threshT < 0)) stop("Threshold must be positive!")
  if(mu<0 || beta<0|| sigma<0) stop("Parameters must be positive!")
  if(nC<1 || round(nC)!=nC) stop("Number of contributors must be a positive whole number!")
  if(stutt<0 || stutt>1) stop("Stutter proportion must be between 0 and 1!")
  if(stuttFW<0 || stuttFW>1) stop("Forrward stutter proportion must be between 0 and 1!")
  if( any(prC<0) || any(prC>1)) stop("Dropin probability must be between 0 and 1!")
  nL<-length(popFreq) #number of loci
  locs <- names(popFreq) #get loci names
  if(is.null(locs)) locs = paste0("locus",1:nL)

  isNotSimplex = function(x) { #helpfunction of consider being simplex
 	 return( round(sum(x),6)!=1 || any(x < 0 | x > 1) )
  }
  
  for(loc in locs) {
    freq = popFreq[[loc]]
    if(isNotSimplex(freq)) warning( paste0(loc," was not a valid simplex") )
	popFreq[[loc]] = freq/sum(freq) #rescaling
  }
  
  #Preparet mixture proportions (mx)
  if(is.null(mx)) {
   mx <- rgamma(nC,1) #generate a flat distribution of mx if not provided
   mx = mx/sum(mx) #rdirichlet(1,rep(1,nC))  #simulate mx for contributors
  } else {
   if( isNotSimplex(mx)) stop("mx is not a simplex")
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

  #convert (mu,sigma) to gamma-parameters
  rho <- 1/(sigma^2)
  tau <- mu/rho

  #for each replicates:
  samples <- list()
  for(r in 1:nrep) { #for each replicate
   mixData <- list()
  # nDropout <- nDropin <- nStutter <- rep(NA,nL) #counts dropout/dropin/stutters (not used)

   for(loc in locs) {  #for each loci
#loc=locs[1]
    #INCLUDE REFERENCE ALLELES OR GENERETATE NEW IF NOT PROVIDED
    if(is.null(refData[[loc]])) {
       refData[[loc]] <- list()     
       nR <- 0 
    } else {
       nR <- sum(sapply(refData[[loc]],function(x) length(x)>0)) #length(refData[[loc]]) #number of non-empty refs.
    }
    mixA <- unlist(refData[[loc]]) #vectorize
    if(nR<nC) { #sample more alleles if missing
     ran <- (nR+1):nC #new range of genref
     for(s in ran) {
      Asim <- refData[[loc]][[s]] <-  sample(names(popFreq[[loc]]),size=2,prob=popFreq[[loc]],replace=TRUE)
      mixA = c(mixA,Asim) #keep not droped
     }
     names(refData[[loc]])[ran] <- paste0("genref",ran)
    }
    mixH = c(t(replicate(2,mx))) #obtain mixture proportions for each alleles
    agg=aggregate(mixH,by=list(mixA),sum) #aggregate contributions for each alleles
    mixA <- agg[,1] #extract allele names

    degscale <- 1 #default is no degradation scale for each allele
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
    
    #Simulate 'TrueAllele' peak heights from gamma distribution:
    mixH <- rgamma(length(mixA),shape=rho*agg$x*degscale,scale=tau) #shape/scale given. We round to integer later.

    #Simulate stutters (only if alleles can be converted to numerics)
    if(stutt>0 || stuttFW>0 ) { #if stutter proportions different fomr zero
     isNum = suppressWarnings( { !any(is.na(as.numeric(mixA))) } ) #check if all alleles can be converted to numbers
     
     if(isNum) {
       mixA2 <- as.numeric(mixA)
       mixA2 <- c(mixA,mixA2-1,mixA2+1) #get stutter positions (original,BW,FW)
       mixH2 <- c( (1-(stutt+stuttFW))*mixH,stutt*mixH,stuttFW*mixH) #get true allele loss, backward/forward stutter contribution
       agg2 <- aggregate(mixH2,by=list(mixA2),sum) #aggregate stutter peaks wrt unique alleles
       mixA <- agg2[,1] #obtain alleles again
       mixH <- agg2$x #obtain peak heights again
     } 
    }
    
    #Get possibly marker specific values threshT,prC,lambda :
    threshT0 = getMarkerVal(threshT,loc) #get relevant detection threshold
    prC0 = getMarkerVal(prC,loc) #get relevant dropin prob
    lam0 = getMarkerVal(lambda,loc) #get relevant dropin rate
    
    #Simulate drop-in
    pos <- ceiling(log(1e-16)/log(prC0)) #get number of possible dropins
    if(pos>1) { #dropin done last
     prCvec <- c(1-prC0/(1-prC0),prC0^(1:pos))
     nDI <- sample(0:pos,size=1,prob=prCvec,replace=TRUE) #number of dropin
     if(nDI>0) {
      DIA <- sample(names(popFreq[[loc]]),size=nDI,prob=popFreq[[loc]],replace=TRUE) #random from population
      if(lam0<=0) stop("Lambda was not specified greater than zero")
      DIH <- rexp(nDI,rate=lam0) + threshT0 #peak height of same alleles are accumulated
      mixH2 <- c(mixH,DIH)
      mixA2 <- c(mixA,DIA)
      agg2=aggregate(mixH2,by=list(mixA2),sum) #aggregate dropped in alleles (assume additative)
      mixA <- agg2[,1]
      mixH <- agg2$x
     }
    }

    dropped <- mixH < threshT0 #alleles to remove because of small PHs (notice locus specific thresholds)
    mixData[[loc]] <- list(adata=mixA[!dropped],hdata=round(mixH[!dropped]))
   } #end each locus
   names(mixData) <- locs
   samples[[r]] <- mixData
  } #end for each rep
  names(samples) <- paste0("sample",1:nrep)
  names(refData) <- locs
  return(list(theta=list(mx=mx,rho=rho,tau=tau,beta=beta,xiBW=stutt,xiFW=stuttFW),samples=samples,refData=refData))
}


