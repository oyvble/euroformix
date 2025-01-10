#' @title plotTopEPG2
#' @author Oyvind Bleka
#' @description EPG data visualizer (interactive)
#' @details Plots the expected peak heights of the top genotypes. The peak heights for corresponding alleles (one sample) are superimposed.
#' @param MLEobj An object returned from calcMLE/contLikMLE
#' @param DCobj An object returned from deconvolve: Must be run with same object as MLEobj
#' @param kit Short name of kit: See supported kits with getKit(). Argument ignored if degradation model used.
#' @param dyeYmax Whether the Y-axis scale should be dye specific
#' @param plotRepsOnly Whether only replicates are plotted 
#' @param options Layout options can be provided: w0 (layout width),marg0 (margin), txtsize0 (text size), locsize0 (locus name size), minY (minimum Y-axis height), ymaxscale (Scaling maximum Y-axis height), transdeg (transparanct of bars 0=max, 1=min)
#' @param withStutterModel Whether taking the stutter model into account when visualizing
#' @return A plotly widget
#' @export

#DCobj=NULL;dyeYmax=TRUE;plotRepsOnly=TRUE;options=NULL;withStutterModel=TRUE
plotTopEPG2 <- function(MLEobj,DCobj=NULL,kit=NULL,dyeYmax=TRUE,plotRepsOnly=TRUE,options=NULL,withStutterModel=TRUE) {
  if(is.null(options$h0)) { h0 = 1200  } else { h0 = options$h0 }  # 5500/nrows #standard height for each dye (depends on number of rows? No)
  if(is.null(options$w0)) { w0 = 1800 } else { w0 = options$w0 }  # standard witdh when printing plot
  if(is.null(options$marg0)) { marg0 = 0.02  } else { marg0 = options$marg0 } 
  if(is.null(options$txtsize0)) { txtsize0 = 15 } else { txtsize0 = options$txtsize0 }  
  if(is.null(options$locsize0)) { locsize0 = 20 } else { locsize0 = options$locsize0 } 
  if(is.null(options$minY)) { minY = 100 } else { minY = options$minY }  #default minimum Y-axis length
  if(is.null(options$ymaxscale)) { ymaxscale = 1.05 } else { ymaxscale = options$ymaxscale }  #default minimum Y-axis length
  if(is.null(options$transdeg)) { transdeg = 0.5 } else { transdeg = options$transdeg }  #transparancy of bars/peaks
  
  #GRAPHICAL SETUP BASED ON SELECTED KIT:
  if(!is.null(MLEobj$model$kit)) kit = MLEobj$model$kit #NB: Use kit of degradation model if given (ignoring kit argument!)
  kitinfo = euroformix::getKit(kit) #names(kitinfo)
  if( is.na(kitinfo)[1] ) {
    print("The kit name was not recognized by getKit!")
    return(NULL)
  }
  kitinfo$Color[kitinfo$Color=="yellow"] = "orange" #exchange col because of better visualization
  
  #Prepare settings:
  dyes = unique(kitinfo$Color)
  nrows = length(dyes) #number of dyes/rows
  
  #Create list with dye,marker,bp (for observed data
  bprng = range(kitinfo$Size) #get range (same range for all plots)
  # bprng[1] = bprng[1]/2 #widen out on left?
  POS_markerNames = aggregate(kitinfo$Size,by=list(kitinfo$Color,kitinfo$Marker),FUN=mean) #get marker positions (bp)
  
  col0 = "#0000007F" #color of all peak heights (is gray colored)
  col1 = .getContrCols(0.3) #get contribution colors
  mx = MLEobj$fit$thetahat2[grep("Mix",names(MLEobj$fit$thetahat2))]
  mxtxt = paste0("Contr.C",1:length(mx),"(",names(col1)[1:length(mx)],")=",signif(mx,3))
  AT = MLEobj$model$AT
  refNames = MLEobj$prepareC$refNamesCond #obtain reference name of conditionals
  nrefs = length(refNames)
  
  #Obtain expected contribution based on model
  df = .getDataToPlotProfile(MLEobj,DCobj,kit,withStutterModel)
  #bprng = range(df$bp)
  sampleNames = unique(df$Sample)
  nReps = length(sampleNames)
  POS_markerNames
  
  #Obtain global maximum window to show
  #Ymax = max(df$PIup)
  Ymax <- max(as.numeric(unlist( strsplit(df$EXP,"/") )))
  Ylim_max =  ymaxscale*max(Ymax,minY,df$Height) #global max y
  
  #SEPARATE PLOTS
  if(nReps==1 || !plotRepsOnly) { #plot separate plot only in this case
  
    for(sample in sampleNames) { #create a seperate EPG plot for each samples
# sample =sampleNames[1]
     plist = list() #create plot object for each color
    
     for(dye in dyes) {
#    dye=dyes[1]
      dyeidx = which(dyes==dye)
      loctab = POS_markerNames[POS_markerNames[,1]==dye,-1,drop=FALSE]
      locs = toupper(as.character(loctab[,1])) #get locs 
      poslocs = loctab[,2] #get corresponding positions
      
      AT1 <- AT #temporary on analytical threshold
      if(!is.null(AT) && length(AT)>1 ) AT1 = AT[ toupper(names(AT))%in%locs ][1]  #extract dye specific AT
      
      #obtain df for specific dye (several markers)
      dfDye = df[df$Sample==sample & df$Marker%in%locs,] #extract subset 
      if(dyeYmax) ymax1 = ymaxscale*max( na.omit( c(minY,AT1,as.numeric(unlist( strsplit(dfDye$EXP,"/") )),dfDye$Height) ) )  #get max 
    
      p = plotly::plot_ly(colors=col0,mode="lines",height=h0) #df,x = ~bp,y=~Height,type="scatter",mode="markers",colors=markerCol,name=~Allele)
      for(j in 1:nrow(dfDye))  p = plotly::add_trace(p,x =dfDye$bp[j] + 1*c( -1/4,0,1/4),y =c(0,dfDye$Height[j],0 ),name=as.character(dfDye$Allele[j]),type = "scatter" , mode = "lines", fill = "tozeroy",fillcolor=col0,showlegend = FALSE,color=factor(1))
      if(!is.null(AT1)) p <- plotly::add_lines(p,x = bprng, y = rep(AT1,2),color=factor(1),line=list(dash = 'dot',width=2),showlegend = FALSE)
      
      for(loc in locs) { #for each loci: color name for different probabilities 
        pGmarker = dfDye$genoProb[dfDye$Marker==loc] #genotype prob (whole marker)
        markerCol = "black" 
        if( length(pGmarker)>0 ) {
          probGeno = min(as.numeric(pGmarker)) 
          if( probGeno>=.95 ) {
            markerCol = "forestgreen"
          } else if(probGeno>=.9) {
            markerCol = "orange"
          } else {
            markerCol = "red"
          }
        }
        p = plotly::add_annotations(p, x=poslocs[which(loc==locs)] ,y=ymax1,text=loc,showarrow=FALSE,font = list(color = markerCol,family = 'Gravitas One',size = locsize0))  #ADD LOCI NAMES
      } #end for each loci
      p = plotly::add_annotations(p, x=dfDye$bp,y=rep(0,nrow(dfDye)),text=dfDye$Allele,showarrow=FALSE,font = list(color = 1,family = 'sans serif',size = txtsize0),yshift=-10)  #ADD ALLELE NAMES
      if(nrefs>0) {
         p = plotly::add_annotations(p, x=dfDye$bp,y=rep(0,nrow(dfDye)),text=dfDye$reftxt,showarrow=FALSE,font = list(color = 1,family = 'sans serif',size = txtsize0),yshift=-25)  #ADD ALLELE NAMES
         if(dyeidx ==1) p = plotly::add_annotations(p, x=rep(bprng[2],2),y=c(ymax1-ymax1/10*(1:nrefs)),text= paste0("Label ",1:nrefs,": ",refNames),showarrow=FALSE,font = list(colors = 1,family = 'sans serif',size = 15),xshift=0,xanchor = 'right')  #ADD ALLELE NAMES
      }
      if(dyeidx ==length(dyes)) p = plotly::add_annotations(p, x=rep(bprng[2],2),y=c(ymax1-ymax1/10*(1:length(mx))),text=mxtxt,showarrow=FALSE,font = list(colors = 1,family = 'sans serif',size = 15),xshift=0,xanchor = 'right')  #ADD ALLELE NAMES
    
      #ADD EXPECTATIONS
      Ev = strsplit(dfDye$EXP,"/") #extract
      shapeList = list()  #add opacity shapes
    
      cc = 1 #counter
      for(i in 1:length(Ev)) {
       Ev2 = c(0,Ev[[i]])
       for(j in 1:length(Ev[[i]])) { #for each contributor 
        if(Ev2[j+1]==Ev2[j]) next #skip if equal
        shapeList[[cc]] = list(type = "rect",fillcolor = col1[j], line = list(color = col1[j],width=0.1), opacity = transdeg,x0 =dfDye$bp[i]-1/3, x1 = dfDye$bp[i]+1/3,y0 = Ev2[j], y1 = Ev2[j+1])#,xref="x", yref = "y")
        cc = cc + 1
       }
      }
      p = plotly::layout(p,xaxis = list(range = bprng,showticklabels = FALSE,title = ""),yaxis=list(range=c(0,ymax1), showline = TRUE,title = "Heights (RFU)"),shapes=shapeList)#,colorway =markerCol) 
      plist[[dye]] = p
     }
     sub = plotly::subplot(plist, nrows = nrows, shareX = FALSE, shareY = FALSE,margin=marg0,titleY= TRUE)
     sub = plotly::layout(sub ,title=sample,barmode = 'group',xaxis = list(title = ""))
     sub = plotly::config(sub, scrollZoom=TRUE, displaylogo=FALSE,modeBarButtonsToRemove=c("hoverClosestCartesian","hoverCompareCartesian","toggleSpikelines"),toImageButtonOptions=list(width=w0)) 
     print(sub) 
    
     if(nReps==1) return(sub) #return function if no replicates
    }
  } #end if not only separate plots
  repcols = rep("gray",nReps)
  
  #REPS IN SAME PLOT
  plist = list() #create plot object for each color
  for(dye in dyes) { #traverse each dye
    #   dye=dyes[1]
    dyeidx = which(dyes==dye)
    loctab = POS_markerNames[POS_markerNames[,1]==dye,-1,drop=FALSE]
    locs = toupper(as.character(loctab[,1])) #get locs 
     
    AT1 <- AT #temporary on analytical threshold
    if(!is.null(AT) && length(AT)>1 ) AT1 = AT[ toupper(names(AT))%in%locs ][1]  #extract dye specific AT
    
    poslocs = loctab[,2] #get corresponding positions
    dfDye = df[df$Marker%in%locs,] #extract df subset 
    dfDye_unique = unique( subset(dfDye,select=c("Marker","Allele","bp","reftxt") ) ) #get unique outcome
    if(dyeYmax) ymax1 = ymaxscale*max(na.omit( c(minY,AT1,dfDye$Height, max(as.numeric(unlist(strsplit(dfDye$EXP,"/")))) ) ))  #get max 
    
    p = plotly::plot_ly(dfDye,type = "bar",height=h0, colors=repcols,showlegend = FALSE)
    p = plotly::add_trace(p,x=~bp,y=~Height,name=~Sample,showlegend = FALSE,color=~Sample,hoverlabel=list(font=list(size=14),namelength=1000,namecolor="black"),text =~Allele) #dfDye,x = ~bp,y=~Height,type="scatter",mode="markers",colors=markerCol,name=~Allele)
    if(!is.null(AT1)) p = plotly::add_segments(p,x = bprng[1], xend = bprng[2], y = AT1, yend = AT1,color=factor(1),line=list(dash = 'dot',width=2),showlegend = FALSE)
    
    for(loc in locs) { #for each loci: color name for different probabilities 
      pGmarker = dfDye$genoProb[dfDye$Marker==loc] #genotype prob (whole marker)
      markerCol = "black" 
      if( length(pGmarker)>0 ) {
        probGeno = min(as.numeric(pGmarker)) 
        if( probGeno>=.95 ) {
          markerCol = "forestgreen"
        } else if(probGeno>=.9) {
          markerCol = "orange"
        } else {
          markerCol = "red"
        }
      }
      p = plotly::add_annotations(p, x=poslocs[which(loc==locs)] ,y=ymax1,text=loc,showarrow=FALSE,font = list(color = markerCol,family = 'Gravitas One',size = locsize0))  #ADD LOCI NAMES
    } #end for each loci
    p = plotly::add_annotations(p, x=dfDye_unique$bp,y=rep(0,nrow(dfDye_unique)),text=dfDye_unique$Allele,showarrow=FALSE,font = list(color = 1,family = 'sans serif',size = txtsize0),yshift=-10)  #ADD ALLELE NAMES
    if(nrefs>0) {
       p = plotly::add_annotations(p, x=dfDye_unique$bp,y=rep(0,nrow(dfDye_unique)),text=dfDye_unique$reftxt,showarrow=FALSE,font = list(color = 1,family = 'sans serif',size = txtsize0),yshift=-25)  #ADD ALLELE NAMES
       if(dyeidx ==1) p = plotly::add_annotations(p, x=rep(bprng[2],2),y=c(ymax1-ymax1/10*(1:nrefs)),text= paste0("Label ",1:nrefs,": ",refNames),showarrow=FALSE,font = list(colors = 1,family = 'sans serif',size = 15),xshift=0,xanchor = 'right')  #ADD ALLELE NAMES
    }
    if(dyeidx ==length(dyes)) p = plotly::add_annotations(p, x=rep(bprng[2],2),y=c(ymax1-ymax1/10*(1:length(mx))),text= mxtxt,showarrow=FALSE,font = list(colors = 1,family = 'sans serif',size = 15),xshift=0,xanchor = 'right')  #ADD ALLELE NAMES
    
    #ADD EXPECTATIONS
    Ev = strsplit(dfDye$EXP,"/") #extract
    shapeList = list()  #add opacity shapes
    cc = 1 #counter
    for(i in 1:length(Ev)) {
     Ev2 = c(0,Ev[[i]])
     for(j in 1:length(Ev[[i]])) { #for each contributor 
      if(Ev2[j+1]==Ev2[j]) next #skip if equal
      shapeList[[cc]] = list(type = "rect",fillcolor = col1[j], line = list(color = col1[j],width=0.1), opacity = transdeg,x0 =dfDye$bp[i]-1/2, x1 = dfDye$bp[i]+1/2,y0 = Ev2[j], y1 = Ev2[j+1])#,xref="x", yref = "y")
      cc = cc + 1
     }
    }
    p = plotly::layout(p,xaxis = list(range = bprng,showticklabels = FALSE,title = ""),yaxis=list(range=c(0,ymax1), showline = TRUE,title = "Heights (RFU)"),shapes=shapeList)#,autosize=FALSE,width=10)#,colorway =markerCol) 
    plist[[dye]] = p
  } #end for each dye
  sub = plotly::subplot(plist, nrows = nrows, shareX = FALSE, shareY = FALSE,margin=marg0,titleY= TRUE) 
  sub = plotly::layout(sub ,title=paste0(sampleNames,collapse="/"),barmode = 'group')
  sub = plotly::config(sub, scrollZoom=TRUE, displaylogo=FALSE,modeBarButtonsToRemove=c("lasso2d","select2d","hoverClosestCartesian","hoverCompareCartesian","toggleSpikelines"),toImageButtonOptions=list(width=w0))
  print(sub)
  
  return(sub)
} #end function
  