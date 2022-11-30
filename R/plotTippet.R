#' @title plotTippet
#' @author Oyvind Bleka
#' @description  Plotting non-contributor plot (tippets)
#' @details This function prompts the user to obtain a filename and saves text based tables and ensures that a proper extension is used (txt or csv)
#' @param x log10 LR values of non-contributors
#' @param lr0 Observed LR to superimpose (log10 LR)
#' @param qq Quantiles to calculate and show in figure
#' @param mtxt Added text to figure
#' @param xrange Range of x-limit (Can be provided by user)
#' @param dig Number of decimals to use
#' @param returnStatsOnly Whether to only return stats
#' @export

plotTippet <- function(x,lr0=NULL,qq=c(0.5,0.95,0.99),mtxt="",xrange=NULL,dig=2,returnStatsOnly=FALSE) {
  
  #Calc summary
  n <- length(x)
  LRavg <- sum(10^x)/n #mean LR
  LRstd <- sqrt( (sum((10^(2*x))) - n*LRavg^2)/(n-1) )  #emperical variance
  quantiles <- quantile(x,qq)
  names(quantiles) = paste0(names(quantiles),"-perc.")
  quantiles <- c(quantiles,Max=max(x))
  
  okLR <- x[!is.infinite(x)]
  n2 <- length(okLR)
  posLR <- n2/n #ratio of positive LRs
  poslogLR <- sum(x>0)/n #ratio of positive logLRs
  
  rateLR=NA #number of observations above observed LR
  if(!is.null(lr0)) rateLR <- sum(x>lr0)/n
  
  #CREATING STAT OBJECT:
  stat = list(posLR,poslogLR,rateLR,quantiles,LRavg,LRstd) #obtain list of statistics
  names(stat[[1]]) = "rate(LR>0)"
  names(stat[[2]]) = "rate(LR>1)"
  names(stat[[3]]) = "rate(LR>obs)"
  names(stat[[5]]) = "Mean LR"
  names(stat[[6]]) = "Std LR"
  if(returnStatsOnly) return(stat)
  
  #Obtain legend text  
  txt = NULL
  for(i in seq_along(stat)) {
    for(j in seq_along(stat[[i]])) {
      var = stat[[i]][j] #obtain variable
      txt = c(txt, paste0( names(var),"=",format(var,digits=dig)))
    }
  }
  
  #print(stat)
  if(n2>5e6) {
    print ("Number of values to plot was above 5mill. This is too large to plot. Look printed values instead.")
  } else if(n2>0) {
    xlim = range(okLR) #obtain range of LR
    if(!is.null(xrange)) {
      xlim[1] <- max(xlim[1],xrange[1]) 
      xlim[2] <- min(xlim[2],xrange[2]) 
    }
    if(!is.null(lr0)) {
      if(lr0>xlim[2]) xlim[2] <- lr0 + 2
      if(lr0<xlim[1]) xlim[1] <- lr0 - 1
    }
    
    plot(ecdf(x), main=paste0("Non-contributor LRs (",n," samples)"),xlab="log10(LR)",xlim=xlim)
    for(q in qq) abline(h=q,col="gray",lty=3)
    abline(v=0,col="gray",lty=2)
    mtext(mtxt)
    if(!is.null(lr0)) { #plot observed LR
      points(lr0,1,pch=10,col="blue")
      lines(rep(lr0,2),c(1,0),lty=1,col="blue",lwd=0.5)
      #print(paste0("Discriminatory metric (log10(LR) - q99) = ",format(lr0 -  quantile(x,0.99),digits=dig) ))
    }
    legend("bottomright",legend=txt,lty=NULL, bg="white",cex=0.8)
  } else {
    print("No positive LR has been simulated")
  } 
  return(stat)
}