#' @title validMCMC
#' @author Oyvind Bleka
#' @description Validates aposteriori samples from MCMC method
#' @details This function takes samples from the MCMC as given in a matrix and shows the aposterior functions.
#' @param mcmcfit A object returned by contLikMCMC
#' @param trace Whether showing trace of samples.
#' @param acf Whether showing autocorrelation function of samples.
#' @param addLastMx Whether showing the mixture proportion of last contributor
#' @export
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
#' dat = prepareData(samples,refData=refData,popFreq=popFreq,threshT=AT) 
#' plotEPG2(dat$samples,dat$refData,kit=kit,AT=AT0)
#' mlefit = contLikMLE(3,dat$samples,dat$popFreq,dat$refData,1:3,
#' 	kit=kit,xi=NULL,prC=0.05,lambda=0.01,seed=1)
#' mcmcfit = contLikMCMC(mlefit,niter=5000,delta=1,seed=1) #mcmcfit$acc
#' validMCMC(mcmcfit,acf=FALSE)
#' }


validMCMC <- function(mcmcfit,trace=TRUE,acf=TRUE, addLastMx = TRUE) {
 mixProptxt = "Mix-prop. C"
 txt <- colnames(mcmcfit$posttheta) #obtain param text
 mxInd = grep(mixProptxt,txt) #obtain mixture propportion indices
 MLEv = mcmcfit$MLE #vector with MLE values
 SDv = sqrt(diag(mcmcfit$Sigma)) #vector with SD values
 if(addLastMx && length(mxInd)>0) { #there must be at least one contributor in order to show
   MLElast =  1-sum(MLEv[mxInd])
   SDlast = sqrt(sum(mcmcfit$Sigma[mxInd,mxInd])) #Var = sum(Cov[mx])
   insMx = max(mxInd) + 1 #Mx inserted last
   
   MLEv = c(MLEv[mxInd],MLElast,MLEv[-mxInd])
   SDv = c(SDv[mxInd],SDlast,SDv[-mxInd])
   txt = c(txt[mxInd],paste0(mixProptxt,insMx),txt[-mxInd])
   
   #insert to posttheta
   posthetaMx = mcmcfit$posttheta[,mxInd,drop=FALSE]
   mcmcfit$posttheta = cbind(posthetaMx, 1-rowSums(posthetaMx),mcmcfit$posttheta[,-mxInd])
 }
 #plot(mcmcfit$postlogL,ty="l")
 
 Ubound <- apply(mcmcfit$posttheta,2,max)
 Lbound <- apply(mcmcfit$posttheta,2,min)
 p <- length(txt)
 par(mfrow=c(p,1+sum(c(trace,acf)) ),mar = c(1.2,1,1,0.2), mgp = c(0,0.2,0))
 for(i in 1:p) {
  dens <- density(mcmcfit$posttheta[,i],from=Lbound[i],to=Ubound[i],n=max(Ubound[i],1024))
  xrange <- range(mcmcfit$posttheta[,i])
  mled <-dnorm(dens$x,MLEv[i],SDv[i] ) #density of lazy bayes
  plot(dens$x,dens$y,ty="l",main=txt[i],xlab="",ylab="",ylim=c(0,max(mled,dens$y)),xlim=xrange )
  lines(dens$x,mled,col=2,lty=2,ylab="",xlab="")
  if(trace) plot(mcmcfit$posttheta[,i],ty="l",ylab="",xlab="")
  if(acf) acf(mcmcfit$posttheta[,i],lag.max=200,ylab="",xlab="")
 }
 dev.new()
 op <- par(no.readonly = TRUE)
 dev.off()
 par(op)
}

