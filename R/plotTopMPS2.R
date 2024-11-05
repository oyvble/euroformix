#' @title plotTopMPS2
#' @author Oyvind Bleka
#' @description MPS data visualizer (interactive)
#' @details Plots the expected peak heights of the top genotypes. The peak heights for corresponding alleles (one sample) are superimposed.
#' @param MLEobj An object returned from calcMLE/contLikMLE
#' @param DCobj An object returned from deconvolve: Must be run with same object as MLEobj
#' @param grpsymbol A separator for each allele giving plot grouping. Useful for separating conventional repeat units (RU) and sequence variant.
#' @param locYmax Whether Y-axis should be same for all markers (FALSE) or not (TRUE this is default)
#' @param options A list of possible plot configurations. See comments below
#' @param withStutterModel Whether taking the stutter model into account when visualizing
#' @return A plotly widget
#' @export

plotTopMPS2 = function(MLEobj,DCobj=NULL,grpsymbol="_",locYmax=TRUE,options=NULL,withStutterModel=TRUE) {
  #MLEobj<<- MLEobj
  if(is.null(options$h0)) { h0 = 300 } else { h0 = options$h0 } # 5500/nrows #standard height for each dye (depends on number of rows? No)
  if(is.null(options$w0)) { w0 = 1800 } else { w0 = options$w0 } # standard witdh when printing plot
  if(is.null(options$marg0)) { marg0 = 0.015 } else { marg0 = options$marg0 } #Margin between subplots
  if(is.null(options$txtsize0)) { txtsize0 = 12 } else { txtsize0 = options$txtsize0 } #text size for alleles
  if(is.null(options$locsize0)) { locsize0 = 20 } else { locsize0 = options$locsize0 } #text size for loci
  if(is.null(options$minY)) { minY = 100 } else { minY = options$minY } #default minimum Y-axis length
  if(is.null(options$ymaxscale)) { ymaxscale = 1.06 } else { ymaxscale = options$ymaxscale } #y-axis scaling to the locus name positions
  if(is.null(options$grptype)) { grptype="group" } else { grptype = options$grptype }#,"stack" "group" is default 
  Qallele = "99"
  
  col0 = "#0000007F" #color of all peak heights (is gray colored)
  col1 = .getContrCols(0.3) #get contribution colors
  mx = MLEobj$fit$thetahat2[grep("Mix",names(MLEobj$fit$thetahat2))]
  mxtxt = paste0("Contr.C",1:length(mx),"(",names(col1)[1:length(mx)],")=",signif(mx,3))
  AT = MLEobj$model$AT
  refNames = MLEobj$prepareC$refNamesCond #obtain reference name of conditionals
  nrefs = length(refNames)
  
  #Obtain expected contribution based on model
  df = .getDataToPlotProfile(MLEobj,DCobj,kit=NULL,withStutterModel,grpsymbol=grpsymbol)
  sampleNames = unique(df$Sample)
  nReps = length(sampleNames)
  locsAll = unique(df$Marker)

  #CREATE PLOTS
  #colnames(df)
  EYmax <- max(as.numeric(unlist( strsplit(df$EXP,"/") )))
  ymax1 <- ymaxscale*max(minY,EYmax,df$Height) #global max y
  
  #GRAPHICAL SETUP BASED ON SELECTED KIT:
  nLocs = length(locsAll) #obtain number of loci
  ncols = 5 #number of locs per row (depend on amax)?
  maxA = max(aggregate(df$Height,by=list(df$Sample,df$Marker),length)$x) #max alleles (any locus)
  if(maxA<=2 && nLocs>40) ncols=10 #Number of cols should depend on allele outcome (SNP vs STR)
  if(!is.null(options$ncols)) ncols=options$ncols #use number of column specified in options 
  h1 = h0*ncols #standard height for each graph (DEPEND ON NUMBER COLS)	
  nrows0=ceiling((nLocs+1)/ncols) #number of rows to use (use 10 per column)
  if(nrefs>0) nrows0=ceiling((nLocs+1)/ncols) #number of rows to use (use 10 per column)
  
  hline <- function(y = 0, xr=0:1) { #helpfunction to draw line
    list(type = "line", x0 = xr[1], x1 = xr[2], y0 = y, y1 = y, line = list(color = col0,dash = 'dot',width=2))
  }
  transdeg = .8 #transparancy degree of bars
  
  #SEPARATE PLOTS
  for(sample in sampleNames) {
#  sample = sampleNames[1]
    plist = list() #create plot object for each marker
    for(loc in locsAll) { #for each locus
#  loc  = locsAll[3]
      AT1 <- AT #temporary on analytical threshold
      if(!is.null(AT) && length(AT)>1 ) AT1 = AT[ toupper(names(AT))==toupper(loc) ]  #extract dye specific AT
      
      dfs = df[df$Sample==sample & df$Marker%in%loc,] #extract subset 
      AlleleCE_unique = unique(dfs$AlleleCE) #get unique alleles
      AlleleCE_nunique = length(AlleleCE_unique) #number of unique
      AlleleCE_unique = sort(AlleleCE_unique,decreasing = FALSE) #sort 
      suppressWarnings({  #Convert CE to numbers (used to sort)
        AlleleCE_uniqueNumeric = as.numeric(AlleleCE_unique)
        if(all(!is.na(AlleleCE_uniqueNumeric))) AlleleCE_unique = sort(AlleleCE_uniqueNumeric,decreasing = FALSE)
      })
      nAlleles = length(dfs$Allele) #number of alleles
      xpos = seq_len(nAlleles)-1 # 0:(nAlleles-1) #position of alleles (standard)
      
      #Obtain number of iso-allels (layers)
      reptab = table(dfs$AlleleCE)
      nLayers = max(reptab) #number of layers (corresponds to number of iso-alleles per CE-allele)
      colsLay = rep(col0,nLayers)  #use defualt color for all layesr
      #repcols = rep("black",nLayers) #gray.colors(nLayers, start = 0.3, end = 0.9) #color level will be adjust regarding #layers
      idxLay = rep(NA,nAlleles)#get layer index
      for(aa in AlleleCE_unique) {
        ind = which(dfs$AlleleCE==aa)
        idxLay[ind] = seq_along(ind)
      }
      xpos1 = rep(NA,nAlleles)
      for(i in seq_len(nAlleles)) xpos1[i] = which(dfs$AlleleCE[i]==AlleleCE_unique)-1
      xpos0 = seq_len(AlleleCE_nunique) - 1
      atxtL = nchar(AlleleCE_unique) #get allele length
      
      p = plotly::plot_ly(dfs,height=h1,showlegend = TRUE)
      p = plotly::add_trace(p,type = "bar", x = xpos1,y=~Height,name=idxLay,hoverinfo="y+text",hoverlabel=list(font=list(size=12),namelength=1000),text =~Allele,color=I(colsLay[idxLay]))
      
      #ADD EXPECTATIONS
      shifts = seq(-0.25,0.25,l=nLayers)
      if(nLayers==1) shifts = 0
      
      Ev = strsplit(dfs$EXP,"/") #extract PHexpectation
      shapeList = list()  #add opacity shapes
      if(!is.null(AT)) shapeList[[length(shapeList) + 1]] = hline(AT1,xr=c(-0.5,AlleleCE_nunique-0.5))
      for(i in seq_along(Ev)) { #for each alleel
        x0 = xpos1[i] #get allele position (centered)
        x1 = x0+shifts[idxLay[i]] #layer decides layer pos
        Ev2 = c(0,Ev[[i]])
        for(j in 1:length(Ev[[i]])) { #for each contributor 
         if(Ev2[j+1]==Ev2[j]) next #skip if equal
          shapeList[[length(shapeList) + 1]] = list(type = "rect",fillcolor = col1[j], line = list(color = col1[j],width=0.1), opacity = transdeg,x0 =x1-1/(2*(nLayers+1)), x1 = x1+1/(2*(nLayers+1)),y0 = Ev2[j], y1 = Ev2[j+1])#,xref="x", yref = "y")
        }
      }
      
      #Add threshold lines to shapes
      if(locYmax)  ymax1 = ymaxscale*max( na.omit( c(minY,AT1,as.numeric(unlist(Ev)),dfs$Height) ))  #get max 
      pGvec = dfs$genoProb
      markerCol = "black" 
      if( length(pGvec)>0 ) {
        prob = min(as.numeric(pGvec)) #get smallest marginal prob
        if( prob>=.95 ) {
          markerCol = "forestgreen"
        } else if(prob>=.9) {
          markerCol = "orange"
        } else {
          markerCol = "red"
        }
      }
      p = plotly::add_annotations(p, x=(AlleleCE_nunique-1)/2 ,y=ymax1,text=loc,showarrow=FALSE,font = list(color = markerCol,family = 'Gravitas One',size = locsize0))  #ADD LOCI NAME

      if(max(atxtL)<=5) p = plotly::add_annotations(p, x=xpos0 ,y=rep(0,AlleleCE_nunique),text=AlleleCE_unique ,showarrow=FALSE,font = list(color = 1,family = 'sans serif',size = txtsize0),yshift=-10)  #ADD ALLELE NAMES
      #Rotate alleles if many alleles? Using tickangle: layout(xaxis=list(tickangle=-45))
      #p = plotly::add_annotations(p, x=(AlleleCE_nunique-1)/2 ,y=ymax1,text=loc,showarrow=FALSE,font = list(color = 1,family = 'Gravitas One',size = locsize0))  #ADD LOCI NAMES 
      
      if(nrefs>0) { #need to take care of different layers (shift on x-axis) + missing PHs
        refuse = which(dfs$reftxt!="")  #!duplicated(dfs$AlleleCE) #only put reference under relevant alleles
        for(rr in refuse) { #for each refs (some are missing PH)
          x0 = xpos1[rr] #get allele position (centered)
          x1 = x0+shifts[idxLay[rr]] #layer decides layer pos
          p = plotly::add_annotations(p, x=x1,y=0,text=dfs$reftxt[rr],showarrow=FALSE,font = list(color = colsLay[idxLay[rr]],family = 'sans serif',size = txtsize0),yshift=-25)  #ADD ALLELE NAMES
          if(dfs$Height[rr]==0)  p = plotly::add_trace(p,type = "scatter",mode="markers", x = x1,y=0,name=idxLay[rr],hoverinfo="y+text",hoverlabel=list(font=list(size=12),namelength=1000),text=dfs$Allele[rr],color=as.factor(idxLay[rr]))  #add a point to missing alleles (to get hovering)
        } #end for each refuse
      } #end if references
      p = plotly::layout(p,xaxis = list(showticklabels = FALSE,title = ""),yaxis=list(range=c(0,ymax1), showline = TRUE,title = ""),shapes=shapeList)
      plist[[loc]] <- p
    } #end for each locus
    # subplot(plist, nrows = nrows0, shareX = FALSE, shareY = FALSE,margin=marg0,titleY= TRUE)%>%plotly::layout(title=sample,barmode = grptype)%>%hide_legend()
    
    #In last plot we will show the Mx and the contributors 
    p = plotly::plot_ly(x = xpos,y=rep(0,length(xpos)),height=h1,showlegend = FALSE,type="scatter",mode="markers")
    p = plotly::add_annotations(p, x=tail(xpos,1),y= c(ymax1-ymax1/6*(1:length(mx))),text= paste0("Contr.C",1:length(mx),"(",names(col1)[1:length(mx)],")=",signif(mx,3)),showarrow=FALSE,font = list(family = 'sans serif',size = 15),xshift=0,xanchor = 'right') 
    
    if(nrefs>0) p = plotly::add_annotations(p, x=0,y=c(ymax1-ymax1/6*(1:nrefs)-ymax1/12),text= paste0("Label ",1:nrefs,": ",refNames),showarrow=FALSE,font = list(family = 'sans serif',size = 15),xshift=0,xanchor = 'left')  #ADD ALLELE NAMES
    p = plotly::layout(p,xaxis = list(showline=FALSE, showticklabels = FALSE,title = ""),yaxis=list(range=c(0,ymax1), showline=FALSE,showticklabels = FALSE,title = ""))#,colorway =dye2) 
    plist[[nLocs+1]] = p
    
    sub = plotly::subplot(plist, nrows = nrows0, shareX = FALSE, shareY = FALSE,margin=marg0,titleY= TRUE)
    sub = plotly::layout(sub, title=sample,barmode = grptype)
    sub = plotly::hide_legend(sub)
    sub = plotly::config(sub, scrollZoom=TRUE, displaylogo=FALSE,modeBarButtonsToRemove=c("lasso2d","select2d","hoverClosestCartesian","hoverCompareCartesian","toggleSpikelines"),toImageButtonOptions=list(width=w0))
    
    print(sub)
  } #end for each samples
  return(sub) #return last created

} #end function
