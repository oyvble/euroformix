#' @title plotSumPH
#' @author Oyvind Bleka
#' @description Plotting sumPH per marker
#' @details The function is fitting gamma distribution to the sumPH per marker to investigate fit: Degradation or sample differences
#' @param samples A List [[sample]][[locus]] with samples which for each samples has locus-list elements with list elements adata and hdata. 'adata' is a qualitative (allele) data vector and 'hdata' is a quantitative (peak heights) data vector.
#' @param kit shortname of kit: Obtained from getKit()
#' @param ignoreMarkers Give string subset of markers that should be ignored in visualization
#' @param LUSsymbol Used to recognize CE in LUS alleles
#' @param MPSsymbol Used to recognize CE in MPS alleles
#' @export


#ignoreMarkers = c("Y","AM");LUSsymbol = "_"; MPSsymbol = ":"
plotSumPH = function(samples,kit=NULL, ignoreMarkers = c("Y","AM"), LUSsymbol = "_", MPSsymbol = ":") {
  sampleNames = unique(names(samples))
  noKit = is.null(kit) || is.na(kit)
  
  ignoreFun = function(loc) any(sapply(c("Y","AM"),function(x) grepl(x,toupper(loc) )))
  locs <- unique(unlist(lapply(samples,function(x) names(x))))
  locs = locs[!locs%in%c("")]
  locUse = locs[!sapply(locs,ignoreFun)] #remove ignored markers

  if(noKit) { #UPDATED: if kit is not recognized: Make simple Sum of peak heights plot
    locs = locUse
    dat <- matrix(0,ncol=length(sampleNames),nrow=length(locs),dimnames = list(locs,sampleNames))
    for(sample in sampleNames) {
      for(loc in locs) {
        if( !loc%in%names(samples[[sample]]) ) next;
        dat[loc,sample] <- sum(samples[[sample]][[loc]]$hdata)  #take sum
      }
    }
    plot(0,0, ylim=c(0,max(dat,na.rm=T)),xlim=c(0,length(locs)),ty="n",main="Summed intensity per loci",axes=F,xlab="",ylab="Intensity")
    axis(2)
    axis(1,at=1:length(locs)-1,labels=locs,las=2,cex.axis=0.7)
    
    mtext(paste0(sampleNames,collapse="/"))
    dw <- 0.15 #width between bars       
    rw <- (1 - dw)/length(sampleNames) #rectange widths
    for(loc in locs) {
      i <- which(loc==locs)
      for(sample in sampleNames) {
        yval <- dat[loc,sample]  #take sum
        if(is.na(yval)) next;
        j <- which(sample==sampleNames)
        low <- i + (j-1)*rw - rw/2
        up <- i + (j-1)*rw + rw/2
        rect(low-1,0,up-1,yval,col=j)
      }
    }
    yd <- c(dat)
    yd <- yd[!is.na(yd)]
    
    checkgamma = function(y,th) {
      xz <- ppoints(length(y))
      qq1 <- qgamma(xz ,shape=2/th[2]^2,scale=th[1]*th[2]^2)
      dev.new() #create new plot
      plot(qq1,sort(y),main="QQ plot of observed intensities (summed)",xlab="Theoretical intensities (gamma distributed)",ylab="Observed intensities")
      abline(a=0,b=1)
      mtext(paste0(sampleNames,collapse="/"))
    }
    tryCatch({ checkgamma(yd,fitgammamodel(yd)) }, error = function(e) e)
    
  } else { #end not kit found
    kitinfo = getKit(kit) #get kitinfo
    
    #Step1: Create df to do degradation-regression with:
    dyes <- unique(kitinfo$Color)
    dat <- NULL
    for(dye in dyes) {
      locs <- toupper(unique(subset(kitinfo$Marker,kitinfo$Color==dye)))
      for(loc in locs) {
        if(!loc%in%locUse) next #skip
        subK <- subset(kitinfo,kitinfo$Color==dye & toupper(kitinfo$Marker)==loc) 
        for(sample in sampleNames) { #for each replicate
          av  <- samples[[sample]][[loc]]$adata #get alleles
          if(is.null(av) || length(av)==0) next
          if(all(grepl(LUSsymbol,av))) { #in case of LUS. Extract only first allele. OK for general
            av  <- sapply(strsplit(av,LUSsymbol),function(x) x[1]) 
          } else if(all(grepl(MPSsymbol,av))) {  #in case of MPS-SEQ. Extract only first allele. OK for general
            av  <- sapply(strsplit(av,MPSsymbol),function(x) x[1])
          }
          avgsize <-  mean(subK$Size[subK$Allele%in%av],na.rm=TRUE) #get average sizedata for alleles
          if(length(avgsize)==0) next
          dat <- rbind(dat, c(sum(samples[[sample]][[loc]]$hdata),avgsize,dye,loc,sample)) #use average for each locus
        } #end for each replicate 
      }  #end for each mixsel
    }
    if(length(dat)==0) {
      gWidgets2::gmessage("The data and kit selected did not match!",title="Incompatible data found",icon="warning")
      return() #INCOMPATIBLE DATA WITH KIT FOUND
    }
    #dat[dat[,3]=="yellow",3] <- "orange"
    #dyes[dyes=="yellow"] <- "orange"
    
    colnames(dat) = c("sumPH","avgSize","Dye","Marker","Sample")
    dat = as.data.frame(dat)
    dat$sumPH = as.numeric(dat$sumPH)
    dat$avgSize = as.numeric(dat$avgSize)
    dat$avgSize2 = (dat$avgSize - 125)/100 #scale fragment length

    srange = range(kitinfo$Size)
    xz <- seq(srange[1],srange[2],l=30)    
    plot(0,0,xlim=srange ,ylim=c(0, max(dat$sumPH)),ty="n",ylab="Sum of the peak heights (rfu) per marker",xlab="Average fragment length",main="Peak height summaries")
    mtext(paste0(sampleNames,collapse="/"))
    for(s in seq_along(sampleNames)) {
      sdat <- dat[dat$Sample==sampleNames[s],,drop=FALSE]
      if(nrow(sdat)==0) next 
      points(sdat$avgSize, sdat$sumPH,col=s,pch=16,cex=1)
      label <- substr(sdat$Marker,0,4)
      text(sdat$avgSize, sdat$sumPH, labels=label,adj=c(0,-0.5),col=s,cex=0.7)
      
      tryCatch({  #avoid crash
        #fit <- lm(log(sdat$sumPH)~sdat$avgSize)
        #lines(xz,exp(fit$coef[1]+xz*fit$coef[2]),col=col,lty=2) 
        lo <- loess(sdat$sumPH~sdat$avgSize)
        lines(xz, predict(lo,xz), ,col=s,lty=2,ty="l",lwd=1.5)
      }, error = function(e) e)
    }
    legend("bottomleft",sampleNames,bty="n",col=seq_along(sampleNames),cex=0.8,pch=19)  
    
    plotquant <- function(th,alpha=0.01) {
      xz <- seq(srange[1],srange[2],l=1000)     
      lines(xz,qgamma(1-alpha/2,shape=2/th[2]^2*th[3]^((xz-125)/100),scale=th[1]*th[2]^2),col="gray")
      lines(xz,qgamma(alpha/2,shape=2/th[2]^2*th[3]^((xz-125)/100),scale=th[1]*th[2]^2),col="gray")
      lines(xz,2*th[1]*th[3]^((xz-125)/100),col="black")
      legend("topright",legend=c("Expectation",paste0(1-alpha,"-coverage")),col=c("black","gray"),lty=1)
    }

    #fit data based on the gamma-model (it may fail):
    tryCatch({ plotquant(th=fitgammamodel(dat$sumPH,dat$avgSize,delta=0.5)) }, error = function(e) {cat("ERROR :",conditionMessage(e), "\n")})
    
    return(NULL) #don't calculate p-values
    tryCatch({
      pvec <- rep(NA,nrow(dat))
      for(i in 1:nrow(dat) ) {
        th <- fitgammamodel(dat$sumPH[-i],dat$avgSize[-i],delta=0.5)
        pvec[i] <- pgamma(dat$sumPH[i],shape=2/th[2]^2*th[3]^(dat$avgSize2[i]),scale=th[1]*th[2]^2) #get probabilities
      }
      pvec[pvec>0.5] <- 1-pvec[pvec>0.5] #make all values smaller than 0.5
      alpha <- 0.05 #signif level
      ind <- which(pvec<(alpha/length(pvec))) #find flagged markers which are below bonferroni-corrected signif
      if(length(ind)>0) {
        text(dat$avgSize[ind],dat$sumPH[ind], labels=format(pvec[ind],digits=2),adj=c(0,1.2),cex=0.7)
      } #end if flaggings
    }, error = function(e) e) #end error catch
  } #end if kit was found

}