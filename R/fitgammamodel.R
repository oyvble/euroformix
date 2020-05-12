#' @title fitgammamodel
#' @author Oyvind Bleka
#' @description fitgammamodel is a function to fit sumPH of the observed peak heights to a gamma model
#' @details
#' This function is used for fitting a gamma regression model
#'
#' @param x a vector of base pairs
#' @param y a vector of sum peak heights per loci
#' @param DEG a boolean of whether degradation (i.e. the bp at x-axis should be considered) 
#' @param niter number of random samples 
#' @param delta Standard deviation of normal distribution when drawing random startpoints. Default is 1.
#' @param plott a boolean of whether plotting fitted model along with data
#' @param alpha The alpha/2 and 1-alpha/2 percentiles will be plotted if x is provided
#' @param offset a number giving shift of bp in degradation model.
#' @param scale a number giving scaling of shifted bp in degradation model.
#' @param th0 Provided start values to be used in optimizer
#' @return Estimated parameters (mean,coefVar,slopePara)
#' @examples
#' th = c(1000,0.8,0.6)
#' n = 100 #number of samples
#' x = seq(10,300,l=n)
#' y = rgamma(n,shape=2/th[2]^2*th[3]^((x-125)/100),scale=th[1]*th[2]^2)
#' fitgammamodel(y,x,DEG=TRUE,niter=10,plott=TRUE,alpha=0.05)
#' fitgammamodel(y,x,DEG=FALSE,niter=10,plott=TRUE,alpha=0.05)
#' @export

fitgammamodel <- function(y,x=NULL,DEG=TRUE,niter=30,delta=1,plott=FALSE,alpha=0.05,offset=125,scale=100,th0=NULL) {
   if(is.null(x)) DEG=FALSE #degradation set to FALSE if x not given
   isOK = !is.na(y) #get non-NA values
   y = y[isOK] #remove NAs
   if(DEG) x = x[isOK] #remove NAs
   
   #FUNCTION TO BE MINIMIZED
   negloglik <- function(th) {
    th[1:2] <- exp(th[1:2]) #log-transform
    if(DEG) {
      th[3] <- 1/(1+exp(-th[3])) #logit-transform (avoid beta>1)
      val <- -sum(dgamma(y,shape=2/th[2]^2*th[3]^((x-offset)/scale),scale=th[1]*th[2]^2,log=TRUE))
    } else {
      val <- -sum(dgamma(y,shape=2/th[2]^2,scale=th[1]*th[2]^2,log=TRUE))
    }
    if(is.infinite(val)) val <- Inf #.Machine$integer.max 
    return(val)
   }

   if(is.null(th0)) {
	   t0 <- log(c(mean(y)/2,sd(y)/mean(y),0.8)) #expected thetas to start with
	   if(!DEG)  t0 <- t0[1:2] #CV used as start value
   } else {
      t0=th0
   }
   cc <- 1 #
   suppressWarnings({
   bestFoo <- nlm(negloglik, t0 ) #store best foo
   while(cc<niter) { #repeat until accepted fitted
    t1 <- rnorm(length(t0),mean=t0,sd=delta) #generate new
    foo <- nlm(negloglik, t1 )
    if(foo$min<bestFoo$min) bestFoo <- foo #if better    
    cc <- cc + 1 #add
   }
  })
  th = exp(bestFoo$est) #get estimated parameters
  if(DEG && length(th)==3 ) {
    th[3] <- 1/(1+exp(-bestFoo$est[3])) #logit-transform (avoid beta>1), IMPORTANT TO AVOID CRASH IN contLikMLE
  }

  if(plott) { #show plot
   par(mfrow=c(1,1))
   if(!is.null(x) && DEG) {
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

