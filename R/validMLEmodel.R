#' @title validMLEmodel
#' @author Oyvind Bleka
#' @description validMLEmodel makes model validation whether the observed peak heights fits the maximum likelihood fitted gamma distribution.
#' @details The cumulative probability of the observed allele peaks are calculated and compared with a uniform distribution.
#' Function calls a 'cumulative procedure' in C++ which uses the package Boost and performs paralellisation (OpenMP).
#'
#' @param mlefit Fitted object using contLikMLE
#' @param kit Shortname of kit used
#' @param plottitle Maintitle text used in the PP-plot
#' @param alpha The significance level used for the envelope test. Default is 0.01
#' @param createplot Boolean of whether plot should be created
#' @param verbose Boolean of whether printing out information about significant alleles
#' @param maxThreads Maximum number of threads to be executed by the parallelization
#' @return retinfo A dataframe with information from the calculated validation model
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
#' dat = prepareData(samples,refData=refData,popFreq=popFreq,threshT=AT) #obtain data to use for analysis
#' plotEPG2(dat$samples,dat$refData,kit=kit,AT=AT0)
#' mlefit = contLikMLE(3,dat$samples,dat$popFreq,dat$refData,1:3,kit=kit,xi=NULL,prC=0.05,lambda=0.01,seed=1)
#' validMLEmodel(mlefit,kit=kit)
#' }

validMLEmodel <- function(mlefit,kit=NULL,plottitle="PP-plot",alpha=0.01,createplot=TRUE,verbose=TRUE,maxThreads=32) {

  #Obtain kitinfo
  if(!is.null(kit)) {
    kitinfo <- getKit(kit) #get kitinfo
    if(length(kitinfo)==1) {
      print("Wrong kit specified. Continue in function.")
    } 
  }
  
  # plottitle="PP-plot";alpha=0.01;createplot=TRUE;verbose=TRUE;maxThreads=32
  pchmax = 26 #number of possible point types
  model <- mlefit$model #take out assumed model with given data
  c <- mlefit$prepareC #Object already stored in mlefit. returned from prepareC function
  
  #Prepare fixed params:
  nM = c$nM #number of markers to evaluate
  ATv = model$threshT
  pCv = model$prC
  lambdav = model$lambda
  fstv = model$fst
  
  #Prepare estimated params:
  nC = model$nC #number of contributors
  thhat = mlefit$fit$thetahat
  mixprop = head(thhat,nC-1)# fitted mix proportions
  mu1 = thhat[nC] #fitted PHexp
  sig1 = thhat[nC+1] #fitted PHvar
  
  #Obtain remaining params
  beta1 = 1 #set default values
  xiB = model$xi #set default values
  xiF = model$xiFW #set default values
  indsel = nC + 2  #get index to insert estimated variable form param vector(thhat)
  if(!is.null(model$kit)) {
    beta1 = thhat[indsel] #get degrad param
    indsel = indsel + 1 #update index
  }
  if(is.null(xiB)) {
    xiB = thhat[indsel] #get BWstutter param
    indsel = indsel + 1 #update index
  }
  if(is.null(xiF)) xiF = thhat[indsel] #get FWstutter param
  
  #Insert to vector
  muv  = rep(mu1,nM)
  sigmav =  as.numeric(rep(sig1,nM))
  betav =  rep(beta1,nM) 
  xiBv = rep(xiB,nM) 
  xiFv = rep(xiF,nM) #always zero
  
  #Obtaining number of alleles/signif level and upper peak height boundary to integrate up to (maxY)
  emptyvec = as.numeric(rep(0,sum(c$nReps*c$nA))) #length of all alleles
  nAtotObs = sum(sapply(model$samples,function(x) sapply(x,function(y) length(y$adata)))) #number of observed alllees
  #Obtain large Y:
  alphaQ <- 0.001 #ensure very far out in quantile (used for estimating probs in gamma-distribution).
  alpha2 <- alphaQ/nAtotObs #"bonferroni outlier"
  suppressWarnings({ #don't show missing allele warning
    maxYobs <- max(sapply(model$samples,function(x) sapply(x,function(y) max(y$hdata)))) #max observation
    maxYexp <- qgamma(1-alpha2,2/sig1^2,scale=mu1*sig1^2) #max observation in theory
    maxY <- ceiling(max(maxYobs,maxYexp)) #get max observed
  }) 
  
  #Calling Special function in C++
  obj <- .C("cumvalgammaC",emptyvec,emptyvec,as.numeric(maxY),c$nC,c$NOK,c$knownGind,as.numeric(mixprop),as.numeric(muv),as.numeric(sigmav),as.numeric(betav),as.numeric(xiBv),as.numeric(xiFv),as.numeric(ATv),as.numeric(pCv),as.numeric(lambdav),as.numeric(fstv),c$nReps,c$nM,c$nA,c$YvecLong,c$FvecLong,c$nTypedLong,c$maTypedLong,c$basepairLong,c$BWvecLong,c$FWvecLong,c$nPS,c$BWPvecLong,c$FWPvecLong,as.integer(maxThreads),c$relGind,c$ibdLong,PACKAGE="euroformix")
  UaPH = obj[[1]] #obtain cumulative vals int_0^y 
  UaMAX = obj[[2]] #obtain cumulative vals int_0^inf 
  #  UaPH/UaMAX
  
  #Structuring values
  locs = names(model$popFreq)
  sn = names(model$samples) #get sample names
  iter = 1 #iterator
  outtab = numeric()
  for(m in 1:c$nM) {
    loc = locs[m]
    tmp = model$popFreq[[loc]] #extract popFreq
    an = names(tmp) #obtain allele names
    
    #obtain dye color for locus
    dye <- 1
    if(!is.null(kit) && length(kitinfo)>1) {
      dye <- unique(kitinfo$Color[toupper(kitinfo$Marker)==loc])
      if(length(dye)==0) dye <- 1 #default as 1 if not found
    }    
    
    for(a in 1:c$nA[m]) { #for each allele
      for(r in 1:c$nReps[a]) { #for each reps
        if(c$YvecLong[iter]>=ATv[m]) {
          Pa =  UaPH[iter]/UaMAX[iter] #)/(1-UaT[iter]) #obtain cumulative probability
          row = c(dye,loc,sn[r],an[a],c$YvecLong[iter],Pa) #insert dye,locus, sample name, allele, height, cumulative probability
          outtab = rbind(outtab, row)
        }
        iter = iter + 1
      } #end for each reps
    } #end fore each alleles
  } #end for each markers
  
  dyevec = outtab[,1] #dye name
  locvec = outtab[,2] #Locus names
  svec = outtab[,3] #sample name
  avec = outtab[,4] #allele names
  yvec = as.numeric(outtab[,5]) #PH names
  cumprobi = as.numeric(outtab[,6]) #Cumulative probs
  
  N <- length(cumprobi) #number of peak heights
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
    outtab = cbind(dyevec,locvec,svec,avec,yvec,cumprobi,pval[ord2]) #print list of outliers
    colnames(outtab) = c("Dye","Marker","Sample","Allele","Height","ProbObserved","pvalue")
    print(outtab[outside,],drop=FALSE) #outside envelop
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
  retinfo = data.frame(Dye=dyevec,Marker=locvec,Sample=svec,Allele=avec,Height=yvec,ProbObs=cumprobi,pvalue=pval[ord2],Significant=outside,stringsAsFactors=FALSE)
  return(retinfo)
  
}

