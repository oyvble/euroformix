#' @title validMLEmodel
#' @author Oyvind Bleka
#' @description validMLEmodel makes model validation whether the observed peak heights fits the maximum likelihood fitted gamma distribution.
#' @details The cumulative probability of the observed allele peaks are calculated and compared with a uniform distribution.
#' Function calls a 'cumulative procedure' in C++ which uses the package Boost and performs paralellisation (OpenMP).
#'
#' @param mlefit Fitted object using calcMLE/contLikMLE
#' @param kit This is no longer used.
#' @param plottitle Maintitle text used in the PP-plot
#' @param alpha The significance level used for the envelope test. Default is 0.01
#' @param createplot Boolean of whether plot should be created
#' @param verbose Boolean of whether printing out information about significant alleles
#' @return retinfo A dataframe with information from the calculated validation model
#' @export

#plottitle="PP-plot";alpha=0.01;createplot=TRUE;verbose=TRUE
validMLEmodel <- function(mlefit,kit=NULL,plottitle="PP-plot",alpha=0.01,createplot=TRUE,verbose=TRUE) {
  cols = c("black","#df462a","#466791","#60bf37","#953ada","#5a51dc","#d49f36","#507f2d","#db37aa","#5b83db","#c76c2d","#552095","#82702d","#55baad","#62aad3","#8c3025","#417d61")
  xlab="Expected probabilities"#: Unif(0,1)"
  ylab="Observed probabilities"#: (Pr(Yj<=yj|Y_{-j}<=y_{-j},Yj>=thresh,model))"
  sz <- 1.5
  
  #c <- mlefit$prepareC #Object already stored in mlefit. returned from prepareC function
  repNames = mlefit$prepareC$repNames
  locNames = mlefit$prepareC$markerNames
  tab = mlefit$MLEvalidTable #copy table from fitted object
  if(is.null(tab)) return(NULL) #couldnt find validation results
  
  N <- nrow(tab) #length(cumprobi) #number of peak heights
  alpha2 <- alpha/N #0.05 significanse level with bonferroni correction
  cumunif <- ((1:N)-0.5)/N #=punif((1:N)-0.5,0,N)
  
  #hist(cumprobi)
  #sum(cumprobi<=0.5)/N
  cumprobi = tab$ProbObs
  ord <- order(cumprobi,decreasing = FALSE)
  ord2 <- match(1:length(ord),ord) #get reverse index
  pval <- 1-pbeta(cumprobi[ord],N*cumunif,N-N*cumunif+1) #one sided p-value
  
  #Must indicate the values below the line (Symmetry)
  ind <- which(cumprobi[ord]<cumunif) #those below the line (two-sided p-value)
  pval[ind] <- pbeta(cumprobi[ord],N*cumunif,N-N*cumunif+1)[ind] #calculate pval in opposite direction
  pval[is.nan(pval)] = 0
  #cumprobi<qbeta(alpha/2,N*cumunif,N-N*cumunif+1) | cumprobi<qbeta(1-alpha/2,N*cumunif,N-N*cumunif+1)
  outside <- pval[ord2] < (alpha2/2) #criterion outside region (divide by 2 to get two-sided)
  tab = cbind(tab,pvalue=pval[ord2],Significant=outside) #update table
  
  if(verbose) { #print info:
    #print points outside the bonferroni-adjusted envolopment
    print(paste0("Total number of peak height observation=",N))
    print(paste0("Significance level=",signif(alpha*100,digits=2),"%"))
    print(paste0("Bonferroni-adjusted significance level=",signif(alpha2*100,digits=2),"%"))
    print("List of observations outside the envelope (with the Bonferroni-level):") #two sided check
    print(tab[tab$Significant,],drop=FALSE)  #print list of outliers (outside envelop)
  }
  
  #plot   
  if(createplot) {
    zones <- matrix(c(1,1,1, 2,5,4, 0,3,0), ncol = 3, byrow = TRUE)
    layout(zones, widths=c(0.3,7,1), heights = c(1,7,.75))
    par(xaxt="n", yaxt="n",bty="n",  mar = c(0,0,0,0))   # for all three titles: 
    plot(x=1,y=1,type="n",ylim=c(-1,1), xlim=c(-1,1))   # fig 1 from the layout
    text(0,0,paste(plottitle), cex=2)  
    plot(x=1,y=1,type="n",ylim=c(-1,1), xlim=c(-1,1))  # fig 2
    text(0,0,paste(ylab), cex=2, srt=90)   
    plot(x=1,y=1,type="n",ylim=c(-1,1), xlim=c(-1,1)) # fig 3
    text(0,0,paste(xlab), cex=2)  
    
    # fig 4: Rigth margin shows distr per replicate
    par(mar = c(1,0,0,0))
    
    plot(0,0,xlim=0:1,ylim=0:1,ty="n")
    for(rep in repNames) {
#rep=repNames[1]
      rind = which(repNames==rep)
      ind = tab$Sample==rep #obtain index of certain replicate
      subDat = tab[ind,]
      locind = match(subDat$Marker,locNames) #obtain loci to plot
      
      points( rep(rind/(length(repNames)+1),sum(ind)), subDat$ProbObs,col=cols[rind],cex=sz,pch=1)#(locind-1)%%pchmax)
    }
    rect(0,0,1,1)
    
    # fig 5, finally, the scatterplot-- needs regular axes, different margins
    par(mar = c(1,2,0,.5), xaxt="s", yaxt="s", bty="n")
    
    plot(0,0,ty="n",xlim=0:1,ylim=c(0,1),cex.axis=sz,asp=1)
    segments(x0=0,y0=0,x1=1,y1=1,lwd=1.2)
    
    #Goodness of fit test  #REMOVED: pval <- ks.test(cumprobi, "punif")$p.value
    #Updated block in 0.6.2: DRAW ENVELOPE LINES FOR ORDER STATISTICS UNDER RANDOMNESS:
    xsq <- seq(0,1,l=1000) #Draw envolope lines
    ysq <- c(alpha,alpha2) #quantiles to consider
    for(qq in ysq) { #for each quanitle
#      qq=ysq[1]
      lines(xsq,qbeta(qq/2,N*xsq,N-N*xsq+1),col=which(ysq==qq),lty=2)
      lines(xsq,qbeta(1-qq/2,N*xsq,N-N*xsq+1),col=which(ysq==qq),lty=2)
    }
    legend("bottomright",legend=paste0("1-Envelope-coverage=",c("",paste0(alpha,"/",N,"=")),signif(ysq,2)),col=1:length(ysq),lty=2,cex=1.3)
    
    #  abline(0,1)
    locind = match(tab$Marker,locNames) #obtain loci to plot
    repind = match(tab$Sample,repNames) #get replicate indices
    points(cumunif,tab$ProbObs[ord],cex=sz,col=cols[repind[ord]],pch=1)#(locind[ord]-1)%%pchmax)
    
    legend("topleft",legend=repNames,pch=19,cex=sz,col=cols[1:length(repNames)])
    dev.new()
    op <- par(no.readonly = TRUE)
    dev.off()
    par(op)
  } #end createplot
  
  #return information
  return(tab)
}