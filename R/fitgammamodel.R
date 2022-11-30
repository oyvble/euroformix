#' @title fitgammamodel
#' @author Oyvind Bleka
#' @description A function to fit sumPH of the observed peak heights to a gamma model 
#' @details This function is used for fitting a gamma regression model
#'
#' @param y Vector with sum peak heights (per loci)
#' @param x Vector with fragment lengths (base pair) (only provided for degradation)
#' @param niter Number of random samples 
#' @param delta Standard deviation of normal distribution when drawing random startpoints. Default is 1.
#' @param plott Whether plotting fitted model along with data
#' @param alpha The alpha/2 and 1-alpha/2 percentiles will be plotted if x is provided
#' @param offset The shift of bp in degradation model.
#' @param scale The scaling of shifted bp in degradation model
#' @param restrictDeg Whether to restrict degradation param to be less than 1
#' @return Estimated parameters (mean,coefVar,slopePara)
#' @export
#' @examples
#'\dontrun{
#' th = c(1000,0.8,0.6)
#' n = 100 #number of samples
#' x = seq(10,300,l=n)
#' y = rgamma(n,shape=2/th[2]^2*th[3]^((x-125)/100),scale=th[1]*th[2]^2)
#' fitgammamodel(y,x,plott=TRUE,alpha=0.05)
#'}

#x=NULL;niter=10;delta=1;plott=FALSE;alpha=0.05;offset=125;scale=100;restrictDeg=FALSE
fitgammamodel <- function(y,x=NULL,niter=10,delta=1,plott=FALSE,alpha=0.05,offset=125,scale=100, restrictDeg=TRUE) {
   DEG = TRUE #degradation used by default
   if(is.null(x) || any(is.na(x))) DEG=FALSE #degradation set to FALSE if x not given or any is NA
   isOK = !is.na(y) #get non-NA values
   y = y[isOK] #remove NAs
   if(DEG) x = x[isOK] #remove NAs
   
   #Impute small Y-vals for zero observations (more robust for low-template profiles)
   minY = min(y[y>0]) #get minimum observed (non-zero)
   y[y==0] = minY/2
    
   #FUNCTION TO BE MINIMIZED
   negloglik <- function(th) {
    th <- exp(th) #assume log-transformed
    if(DEG) {
      #th[3] <- 1/(1+exp(-th[3])) #logit-transform (avoid beta>1)
      val <- -sum(dgamma(y,shape=2/th[2]^2*th[3]^((x-offset)/scale),scale=th[1]*th[2]^2,log=TRUE))
    } else {
      val <- -sum(dgamma(y,shape=2/th[2]^2,scale=th[1]*th[2]^2,log=TRUE))
    }
    if( is.infinite(val)) val <- Inf #.Machine$integer.max 
    return(val)
   }
   
   #Helpfunction for optmization (may throw error)
   helpOptimize = function(th) {
      opt = list(min=Inf)
      suppressWarnings({
         tryCatch({ opt = nlm(negloglik, th )}
         , error = function(e) e )
      })
      return(opt)
   }
   
   #Obtain good start values
   mu0 = mean(y/2)
   omega0 = sd(y/2)/mu0 #0.4 #default start value for PHvar  
   if(DEG) { #if degradation 
     coefs = exp(coef(lm(log(y)~x))) #prefit using logged values (convert params back)
     mu0 = coefs[1]/2 #heterozugout allele
     deg0 = coefs[2] #degrad slope
     sigma0 = sqrt(mean((2*mu0*deg0^x - y)^2))
     omega0 = sigma0/mu0
     #xv = seq(min(x),max(x),l=1000)
     #plot(x,y)
     #lines(xv,2*mu0*deg0^xv,ty="l")
     th0 <- c(mu0,omega0,deg0)
   } else { #if no degradation
     th0 <- c(mu0,omega0)
   }
   #if(verbose) cat(paste0("theta0=",paste0(th0,collapse=",")))  
   
   #Provide repeated optimizations
   largeVal = 1e300 #a large value (below inf)
   t0 = log(th0) #log transform pre-estimates
   bestFoo = helpOptimize(t0)
   if(niter==1 && bestFoo$minimum>largeVal) niter=30 #if 1st optimization failed
     cc <- 1 #
     while(cc<niter) { #repeat until accepted fitted
      t1 <- rnorm(length(t0),mean=t0,sd=delta) #generate new
      foo <- helpOptimize(t1)
      if(foo$min<bestFoo$min) bestFoo <- foo #if better    
      cc <- cc + 1 #add
     }
   if(bestFoo$minimum>largeVal) return(NULL) #return NULL if optimim failed
   
  th = exp(bestFoo$est) #get estimated parameters
  if(restrictDeg && DEG && th[3]>1 ) th[3] = 0.999 #set start value of degrad close to 1 (IMPORTANT TO AVOID CRASH IN contLikMLE)

  if(plott) { #show plot
   par(mfrow=c(1,1))
   if(DEG) {
    plot(x,y,main="Gamma regression")
    xv = seq(min(x),max(x),l=100)
    abline(v=offset,lty=2) 
    lines(xv,2*th[1]*th[3]^((xv-offset)/scale),col="black",lwd=2)
    lines(xv,qgamma(1-alpha/2,shape=2/th[2]^2*th[3]^((xv-offset)/scale),scale=th[1]*th[2]^2),col="gray",lwd=2)
    lines(xv,qgamma(alpha/2,shape=2/th[2]^2*th[3]^((xv-offset)/scale),scale=th[1]*th[2]^2),col="gray",lwd=2)
    legend("topright",legend=c("Expectation",paste0("Lower/upper ",alpha/2*100,"%-percentile")),col=c("black","gray"),lwd=2)
   } else { #show distribution only
    h <- hist(y,breaks=25,probability=T)
    lines(h$mids,dgamma(h$mids,shape=2/th[2]^2,scale=th[1]*th[2]^2),col="black",lwd=2)
   }
  }
  return(th)
} 

