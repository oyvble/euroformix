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
#' @return Estimated parameters (mean,coefVar,slopePara)
#' @examples
#' th = c(1000,0.8,0.6)
#' n = 100 #number of samples
#' x = seq(10,300,l=n)
#' y = rgamma(n,shape=2/th[2]^2*th[3]^((x-125)/100),scale=th[1]*th[2]^2)
#' fitgammamodel(y,x,DEG=TRUE,niter=10,plott=TRUE,alpha=0.05)
#' fitgammamodel(y,x,DEG=FALSE,niter=10,plott=TRUE,alpha=0.05)
#' @export

fitgammamodel <- function(y,x=NULL,DEG=TRUE,niter=30,delta=1,plott=FALSE,alpha=0.05,offset=125,scale=100) {
   if(is.null(x)) DEG=FALSE #degradation set to FALSE if x not given
   negloglik <- function(th) {
    th <- exp(th)
    if(DEG) val <- -sum(dgamma(y,shape=2/th[2]^2*th[3]^((x-offset)/scale),scale=th[1]*th[2]^2,log=TRUE))
    if(!DEG) val <- -sum(dgamma(y,shape=2/th[2]^2,scale=th[1]*th[2]^2,log=TRUE))
    if(is.infinite(val)) val <- .Machine$integer.max 
    return(val)
   }
#   t0 <- log(c(mean(y)/2,0.4,0.8)) #expected thetas to start with
#   if(!DEG)  t0 <- c( t0[1],log(sd(y)/mean(y)) ) #CV used as start value
   t0 <- log(c(mean(y)/2,sd(y)/mean(y),0.8)) #expected thetas to start with
   if(!DEG)  t0 <- t0[1:2] #CV used as start value
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

  if(plott) { #show plot
   par(mfrow=c(1,1))
   if(!is.null(x) & DEG) {
    plot(x,y,main="Gamma regression")
    abline(v=125,lty=2) 
    lines(x,2*th[1]*th[3]^((x-125)/100),col="black",lwd=2)
    lines(x,qgamma(1-alpha/2,shape=2/th[2]^2*th[3]^((x-125)/100),scale=th[1]*th[2]^2),col="gray",lwd=2)
    lines(x,qgamma(alpha/2,shape=2/th[2]^2*th[3]^((x-125)/100),scale=th[1]*th[2]^2),col="gray",lwd=2)
    legend("topright",legend=c("Expectation",paste0("Lower/upper ",alpha/2*100,"%-percentile")),col=c("black","gray"),lwd=2)
   } else { #show distribution only
    h <- hist(y,breaks=25,probability=T)
    lines(h$mids,dgamma(h$mids,shape=2/th[2]^2,scale=th[1]*th[2]^2),col="black",lwd=2)
   }
  }
  return(th)
} 

