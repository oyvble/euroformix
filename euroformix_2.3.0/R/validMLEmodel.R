#' @title validMLEmodel
#' @author Oyvind Bleka
#' @description validMLEmodel makes model validation whether the observed peak heights fits the maximum likelihood fitted gamma distribution.
#' @details The cumulative probability of the observed allele peaks are calculated and compared with a uniform distribution.
#' Function calls procedure in C++ by using the package Armadillo and Boost.
#'
#' @param mlefit Fitted object using contLikMLE
#' @param kit Shortname of kit used
#' @param plottitle Maintitle text used in the PP-plot
#' @param alpha The significance level used for the envelope test. Default is 0.01
#' @param createplot Boolean of whether plot should be created
#' @param verbose Boolean of whether printing out information about significant alleles
#' @return retinfo A dataframe with information from the calculated validation model
#' @export

validMLEmodel <- function(mlefit,kit=NULL,plottitle="PP-plot",alpha=0.01,createplot=TRUE,verbose=TRUE) {
 require(cubature) #Necessary to calculate the integrals
 pchmax = 26 #number of possible point types
 theta2 <- theta <- mlefit$fit$thetahat #condition on mle parameter
 np <- length(theta) #number of unknown parameters
 model <- mlefit$model #take out assumed model with given data
 nC <- model$nC #number of contributors
 nodeg <- is.null(model$kit) #check for degradation
 if(nodeg) theta <- c(theta[1:(nC+1)],1) #insert beta variable equal 1
 if(!is.null(model$xi)) {
  theta <- c(theta,as.numeric(model$xi))
 } else {
  if(nodeg)  theta <- c(theta,theta2[np]) #insert fitted xi last again
 }
 locs <- names(model$popFreq)
 nL <- length(locs)
 mvec <- mlefit$fit$thetahat2[1:nC]
 mu <- theta[nC]
 sigma <- theta[nC+1]
 xi <- model$xi 
 fst <- model$fst

if(is.null(xi)) { #stutter is unknown
  likYtheta <- function(yval) {   #call c++- function: length(theta)=nC+2
    ret$obsY[j] <- yval
    Cval  <- .C("loglikgammaC",as.numeric(0),as.numeric(theta),as.integer(np),ret$nC,ret$nK,ret$nL,ret$nS,ret$nA,ret$obsY,ret$obsA,ret$CnA,ret$allASind,ret$nAall,ret$CnAall,ret$Gvec,ret$nG,ret$CnG,ret$CnG2,ret$pG,ret$pA, as.numeric(model$prC), ret$condRef,as.numeric(model$threshT),as.numeric(model$fst),ret$mkvec,ret$nkval,as.numeric(model$lambda),as.numeric(ret$bp),as.integer(0),PACKAGE="euroformix")[[1]]
    return(exp(Cval + log(model$pXi(theta[ret$nC+3])))) #weight with prior of tau and 
  }
} else { #stutter is known. call c++- function: length(theta)=nC+1
  likYtheta <- function(yval) {   
    ret$obsY[j] <- yval
    Cval  <- .C("loglikgammaC",as.numeric(0),as.numeric(theta),as.integer(np),ret$nC,ret$nK,ret$nL,ret$nS,ret$nA,ret$obsY,ret$obsA,ret$CnA,ret$allASind,ret$nAall,ret$CnAall,ret$Gvec,ret$nG,ret$CnG,ret$CnG2,ret$pG,ret$pA, as.numeric(model$prC), ret$condRef,as.numeric(model$threshT),as.numeric(model$fst),ret$mkvec,ret$nkval,as.numeric(model$lambda),as.numeric(ret$bp),as.integer(0),PACKAGE="euroformix")[[1]]
    return(exp(Cval))
  }
}

  delta <- 0.01 #relative error in integral
  alphaQ <- 0.001 #ensure very far out in quantile (used for estimating probs in gamma-distribution).
  alpha2 <- alphaQ/sum(sapply(model$samples,function(x) sapply(x,function(y) length(y$adata)))) #"bonferroni outlier"
  
 suppressWarnings({ #don't show missing allele warning
   maxYobs <- max(sapply(model$samples,function(x) sapply(x,function(y) max(y$hdata)))) #max observation
   maxYexp <- qgamma(1-alpha2,2/sigma^2,scale=mu*sigma^2) #max observation in theory
   maxY <- ceiling(max(maxYobs,maxYexp)) #get max observed
 }) 
  minY <- model$threshT
 
  if(!is.null(kit)) {
   kitinfo <- getKit(kit) #get kitinfo
   if(length(kitinfo)==1)  print("Wrong kit specified.")
  }

  cumprobi <- avec <- numeric()
  locvec <- dyevec <- numeric()
  for(loc in locs) { #traverse for each locus
   samples <- lapply(model$samples,function(x) x[loc])
   ret <- prepareC(nC=model$nC,samples,popFreq=model$popFreq[loc],refData=model$refData[loc],condOrder=model$condOrder,knownRef=model$knownRef,kit=model$kit,model$knownRel,model$ibd,fst,incS=is.null(xi) || xi>0)
   Yupper <- ret$obsY  #observed peak heights is upper limit in integral
   n <- length(Yupper) #number of observed peak heights 
   if(n==0) next #skip if none observed
   avec <- c(avec,unlist(lapply(samples,function(x) x[[loc]]$adata)) )
   for(j in 1:n) {
    ret$obsY <- Yupper #reset observations
    num <- adaptIntegrate(Vectorize(likYtheta),minY,Yupper[j],tol=delta)[[1]]
    denom <- adaptIntegrate(Vectorize(likYtheta),minY,maxY,tol=delta)[[1]]
    val <- num/denom
    cumprobi <- c(cumprobi,val) #get cumulative probability
    locvec <- c(locvec,loc)
    dye <- 1
    if(!is.null(kit) && length(kitinfo)>1) {
     dye <- unique(kitinfo$Color[toupper(kitinfo$Marker)==loc])
     if(length(dye)==0) dye <- 1 #default as 1 if not found
    }    
    dyevec <- c(dyevec,dye) 
   }
  }

  N <- length(cumprobi) #number of samples
  alpha2 <- alpha/N #0.05 significanse level with bonferroni correction
  
  cumunif <- ((1:N)-0.5)/N #=punif((1:N)-0.5,0,N)
  locind <- sapply(locvec,function(x) which(x==locs))
  #sum(cumprobi<=0.5)/N
  ord <- order(cumprobi)
  ord2 <- match(1:length(ord),ord) #get reverse index
  #plottitle <- "PP-plot between fitted model and theoretical model"
  xlab="Expected probabilities"#: Unif(0,1)"
  ylab="Observed probabilities"#: (Pr(Yj<=yj|Y_{-j}<=y_{-j},Yj>=thresh,model))"
  sz <- 1.5
  dyevec[dyevec=="yellow"] <- "orange" 
 
  pval <- 1-pbeta(cumprobi[ord],N*cumunif,N-N*cumunif+1) #one sided p-value
  ind <- cumprobi[ord]<cumunif #those below the line (two-sided p-value)
  pval[ind] <- pbeta(cumprobi[ord],N*cumunif,N-N*cumunif+1)[ind]
 #cumprobi<qbeta(alpha/2,N*cumunif,N-N*cumunif+1) | cumprobi<qbeta(1-alpha/2,N*cumunif,N-N*cumunif+1)
  outside <- pval[ord2] < (alpha2/2) #criterion outside region (divide by 2 to get two-sided)

 if(verbose) { #print info:
  #print points outside the bonferroni-adjusted envolopment
  print(paste0("Total number of peak height observation=",N))
  print(paste0("Significance level=",signif(alpha*100,digits=2),"%"))
  print(paste0("Bonferroni-adjusted significance level=",signif(alpha2*100,digits=2),"%"))
  print("List of observations outside the envelope (with the Bonferroni-level):") #two sided check
  outtab = cbind(dyevec,locvec,avec,cumprobi,pval[ord2])
  colnames(outtab) = c("Dye","Marker","Allele","ProbObserved","pvalue")
  print(outtab[outside,]) #outside envelop
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
  # fig 4:
  par(mar = c(1,0,0,0))
  cols <- unique(dyevec)
  plot(0,0,xlim=0:1,ylim=0:1,ty="n")
  for(i in 1:length(cols)) {
   inds <- cols[i]==dyevec
   points(rep(i/(length(cols)+1),sum(inds)), cumprobi[inds],pch=(locind[inds]-1)%%pchmax,col=cols[i],cex=sz)
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
  for(qq in ysq) {
    lines(xsq,qbeta(qq/2,N*xsq,N-N*xsq+1),col=which(ysq==qq),lty=2)
    lines(xsq,qbeta(1-qq/2,N*xsq,N-N*xsq+1),col=which(ysq==qq),lty=2)
  }
  legend("bottomright",legend=paste0("1-Envelope-coverage=",c("",paste0(alpha,"/",N,"=")),signif(ysq,2)),col=1:length(ysq),lty=2,cex=1.3)
   
#  abline(0,1)
  points(cumunif,cumprobi[ord],pch=(locind[ord]-1)%%pchmax,cex=sz,col=dyevec[ord])

  locinfo <- unique(cbind(locvec,dyevec,locind)) #info to show in legend
  legend("topleft",legend=locinfo[,1],pch=(as.numeric(locinfo[,3])-1)%%pchmax,cex=sz,col=locinfo[,2])
  dev.new()
  op <- par(no.readonly = TRUE)
  dev.off()
  par(op)
 } #end createplot

  #return information
  retinfo = data.frame(Dye=dyevec,Marker=locvec,Allele=avec,ProbObs=cumprobi,pvalue=pval[ord2],Significant=outside,stringsAsFactors=FALSE)
  return(retinfo)
}



