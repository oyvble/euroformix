#' @title validMCMC
#' @author Oyvind Bleka
#' @description Validates aposteriori samples from MCMC method
#' @details This function takes samples from the MCMC as given in a matrix and shows the aposterior functions.
#' @param mcmcfit A object returned by contLikMCMC
#' @param trace Whether showing trace of samples.
#' @param acf Whether showing autocorrelation function of samples.
#' @export

#trace=TRUE;acf=TRUE
validMCMC <- function(mcmcfit,trace=TRUE,acf=TRUE) {
  MLEv = mcmcfit$MLE #vector with MLE values
  SDv = sqrt(diag(mcmcfit$Sigma)) #vector with SD values
  useVars = SDv>0
  #plot(mcmcfit$postlogL,ty="l")
  p <- sum(useVars) #number of parameters to show
  par(mfrow=c(p,1+sum(c(trace,acf)) ),mar = c(1.2,1,1,0.2), mgp = c(0,0.2,0))
  for(i in seq_along(useVars)) {
    if(!useVars[i]) next
    txt = names(MLEv)[i]
    postTheta = mcmcfit$posttheta[,i]
    xrange <- range(postTheta)
    dens <- density(postTheta,from=xrange[1],to=xrange[2],n=max(xrange[2],1024))
    mled <-dnorm(dens$x,MLEv[i],SDv[i] ) #density of lazy bayes
    if(any(is.infinite(mled))) next
    plot(dens$x,dens$y,ty="l",main=txt,xlab="",ylab="",ylim=c(0,max(mled,dens$y)),xlim=xrange )
    lines(dens$x,mled,col=2,lty=2,ylab="",xlab="")
    if(trace) plot(postTheta,ty="l",ylab="",xlab="")
    if(acf) acf(postTheta,lag.max=200,ylab="",xlab="")
  }
  #RESET parmfrow settings
  dev.new()
  op <- par(no.readonly = TRUE)
  dev.off()
  par(op)
}

