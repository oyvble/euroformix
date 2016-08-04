#' @title plotEPG
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description EPG plotter created by Oskar Hansson.
#' @details Plots peak height with corresponding allele for one sample for a given kit.
#' @param Data List of adata- and hdata-elements.
#' @param kit name of kit: {"ESX17","ESI17","ESI17Fast","ESX17Fast","Y23","Identifiler","NGM","ESSPlex","ESSplexSE","NGMSElect","SGMPlus","ESX16", "Fusion","GlobalFiler"}
#' @param sname Sample name label.
#' @param threshT The detection threshold can be shown in gray in the plot.
#' @param refcond condition on a list$refname$locname$adata of reference alleles which are labeled in EPG
#' @export
plotEPG <- function(Data,kitname,sname="",threshT=0,refcond=NULL) {
 #Data is list with allele and height data. Only one sample!
 #for selected sample:

 generateEPG<-function(typingKit, alleleList, peakHeightList, locusVector, sampleName, refLabelList=NULL, refnames=NULL, dyeVector=NULL,  offsetVector=NULL,repeatUnitVector=NULL,  drawBoxPlots=TRUE, drawPeaks=TRUE){

	# This function generates an electropherogram (EPG) like plot from lists of allele names and peak heights.
	# Default values are used if no typing kit is specified.
	# If a typing kit is specified some or all of its properties can be overriden.

	# INPUT PARAMETERS:

	# alleleList - a list of alleles names.
	# peakHeightList - a list of peak heights.
	# dyeVector - a vector specifying the color of markers.
	# locusVector - a vector specifying the marker names.
	# offsetVector - a vector specifying the marker start offsets in base pairs.
	# repeatUnitVector - a vector specifying the repeat unit size in base pairs.
	# sampleName - title on the EPG.
	# refLabelList - a list of reference-index
	# refnames - a vector with reference names
 	# typingKit - kit short name or index number as specified in the function 'getKit'.

	# NB! If 'locusVector', 'dyeVector', 'offsetVector' or 'repeatUnitVector' is provided
	# they override the information in a matching 'typingKit' IF they are of the same lenght.


	# FLAGS:

	kitFound <- FALSE


	# CONSTANTS:

	# Default sample name / epg header:
	defSampleName <- "EPG"

	# Default repeat unit in base pairs:
	defaultRepeatUnit <- 4

	# Default locus name (will be suffixed with a number):
	defaultLocusName <- "Locus "

	# Default marker spacing in base pairs:
	defaultMarkerSpacing <- 100

	# GRAPH CONSTANTS:
	
	# Graph title:
	# sampleName is used as title.
	
	# X axis label:
#	xlabel <- "base pair"
	xlabel <- "bp"

	# Y axis label:
	ylabel <- "peak height (rfu)"	

	# Peak half width in base pair (width of peak = 2 * peakHalfWidth):
	peakHalfWidth <- 0.6

	# Relative Allele name text size:	
	alleleNameTxtSize <- 0.8

	# Distance between the highest peak in a plot and the plot border (1.04 = 4% margin).
	yMarginTop <- 1.04


	# VERIFICATIONS
	# REQUIRED PARAMETERS:	

	# Check if allele name list is provided.
	if(is.null(alleleList) || !is.list(alleleList))	{

		stop("'alleleList' must be a list giving the alleles per marker/locus")
	}

	# Check if peak height list is provided.
	if(is.null(peakHeightList) || !is.list(peakHeightList)) {

		stop("'peakHeightList' must be a list giving the peak heights per marker/locus")
	}

	# VERIFICATIONS
	# OPTIONAL PARAMETERS:	

	# Check if sample name is provided.
	if (is.null(sampleName)) {
		sampleName <- defSampleName
	}

	# Check if kit name is provided.
	if (!is.null(typingKit)){

		# Get kit information.
		kit <- getKit(typingKit)
	
		# Check if kit was found.	
		if (is.na(kit)[1]) {
			print("Kit not found.")
			print("Using default values.")
			kitFound <- FALSE
			
		} else { # Kit found. Save information in vectors.
		
			kitFound <- TRUE

                  unittab <- getKit(typingKit, what="offset")
                  rangetab <- getKit(typingKit, what="range") 
#cbind(dyetab,unittab) 

			# Marker/locus names.
			locusVectorKit <- rangetab[,1]

			# Dye for each marker/locus
			dyeVectorKit <- rangetab[,2] 
                  #dyeVectorKit <- #toupper(substr(dyeVectorKit, 1, 1)) #Convert!

			#Base pair start offset for each marker.
			offsetVectorKit <- unittab[,2]

			#Size in base pair of repeating unit for each marker.
			repeatUnitVectorKit <- unittab[,3]

                  #Range of markers:
                  locusMinVectorKit <- rangetab[,3]
                  locusMaxVectorKit <- rangetab[,4]

		}
	}

	# If kit is specified...
	if (kitFound){
		# ... check length of allele list.
		numberOfAlleles <- length(alleleList)
		numberOfAllelesKit <- length(locusVectorKit)
		if (numberOfAlleles != numberOfAllelesKit) {
			numberOfAlleles <- numberOfAlleles + 1
			# Extend with missing alleles (add NAs)
			for (index in numberOfAlleles:numberOfAllelesKit) {
				alleleList[[index]] <- c(NA)
			}
		}
		# ... and length of peak height list.
		numberOfPeaks <- length(peakHeightList)
		numberOfAllelesKit <- length(locusVectorKit)
		if (numberOfPeaks != numberOfAllelesKit) {
			numberOfPeaks <- numberOfPeaks + 1
			# Extend with missing peaks (add NAs)
			for (index in numberOfPeaks:numberOfAllelesKit) {
				peakHeightList[[index]] <- c(NA)
			}
		}
	}

	# Check dye vector.
	if (is.null(dyeVector)) {
		if (kitFound) {
			dyeVector <- dyeVectorKit
		} else {
			# No dye vector is specified. Use default color.
			dyeVector <- rep("black", length(alleleList))
		}
	} else {
		if (kitFound) {
			# Use kit specifications if incorrect length.
			if (length(dyeVector) != length(dyeVectorKit)) {
				print("'dyeVector' is not of the same length as kit specifications.")
				print("Kit specifications will be used.")
				dyeVector <- dyeVectorKit
			} else {
				# Override kit specifications with provided values.
			}
		}
	}

	# Check locus vector.
	if (is.null(locusVector)) {
		if (kitFound) {
			locusVector <- locusVectorKit
		} else {
			# No locus vector is specified. Use default name + number suffix.
			locusVector <- paste(rep(defaultLocusName, length(alleleList)), seq(1:length(alleleList)))
		}
	} else {
		if (kitFound) {
			# Use kit specifications if incorrect length.
			if (length(locusVector) != length(locusVectorKit)) {
				print("'locusVector' is not of the same length as kit specifications.")
				print("Kit specifications will be used.")
				locusVector <- locusVectorKit
			} else {
				# Override kit specifications with provided values.
			}
		}
	}
	# Check offset vector.
	if (is.null(offsetVector)) {
		if (kitFound) {
			offsetVector <- offsetVectorKit
		} else {
			# No offset vector is specified. Use default marker spacing.
			offsetVector <- seq(0, length(alleleList) * defaultMarkerSpacing, by = defaultMarkerSpacing)
		}
	} else {
		if (kitFound) {
			# Use kit specifications if incorrect length.
			if (length(offsetVector) != length(offsetVectorKit)) {
				print("'offsetVector' is not of the same length as kit specifications.")
				print("Kit specifications will be used.")
				offsetVector <- offsetVectorKit
			} else {
				# Override kit specifications with provided values.
			}
		}
	}
	# Check repeatUnit vector.
	if (is.null(repeatUnitVector)) {
		if (kitFound) {
			repeatUnitVector <- repeatUnitVectorKit
		} else {
			# No repeatUnit vector is specified. Use default repeat unit length.
			repeatUnitVector <- rep(defaultRepeatUnit, length(alleleList))
		}
	} else {
		if (kitFound) {
			# Use kit specifications if incorrect length.
			if (length(repeatUnitVector) != length(repeatUnitVector)) {
				print("'repeatUnitVector' is not of the same length as kit specifications.")
				print("Kit specifications will be used.")
				repeatUnitVector <- repeatUnitVectorKit
			} else {
				# Override kit specifications with provided values.
			}
		}
	}


	# PREPARE PARAMETERS

	# Convert dye vector to numeric vector.
	colorVector <- dyeVector
	colors <- unique(colorVector)  #added ØB
        noColors <-  length(colors) #added ØB

	# INITIATE VARIABLES

	# Create lists:
	basePairList <- list()
	phListByColor <- list()
	allelesByColorList <- list()
	
      markerByColorList <- list()
      bpmarkerByColorList <- list()
      refLabelByColorList <- list()
      isLab <- length(refnames)>0 #boolean have labels

	# SORT DATA ACCORDING TO COLOR CHANNEL and
	# CONVERT ALLELE NAMES TO FRAGMENT LENGTH IN BASE PAIRS

	# Loop over all color channels.
	for (color in 1:noColors){
            basePairTmpLst <- list()
		# Boolean vector indicating selected markers (same color).
 		selectedMarkers <- colorVector == colors[color]

            if (kitFound) {
             markerByColorList[[color]] <- locusVectorKit[selectedMarkers] #EDIT by ØB: marker names
		} else {
             markerByColorList[[color]] <- locusVector #EDIT by ØB: marker names
            }

		# Extract all alleles in the same color channel.
		allelesByColor <- alleleList[selectedMarkers]

		# Extract all peak heights in the same color channel.
		peakHeightsByColor <- peakHeightList[selectedMarkers]

		# Extract all marker offsets in the same color channel.
		offsetByColor <- offsetVector[selectedMarkers]

		# Extract all repeat unit sizes in the same color channel.
		repeatUnitByColor <- repeatUnitVector[selectedMarkers]

            refLabelByColor <- refLabelList[selectedMarkers] #added ØB

		# Loop over all markers in that color channel.
		for (marker in 1:length(allelesByColor)){

			# Get alleles for the current marker.
			alleleValue <- allelesByColor[[marker]]
			# Check presence of X/Y.
			indexOfXY <- grep("[X,x,Y,y]", alleleValue)
			if (length(indexOfXY)) {
				alleleValue <- toupper(alleleValue)
				# Use 1 and 2 for X and Y.
				alleleValue <- sub("X", 1, alleleValue)
				alleleValue <- sub("Y", 2, alleleValue)
				alleleValue <- as.numeric(alleleValue)
			} else { #added ØB
			 alleleValue <- as.numeric(alleleValue) #added ØB
                  } #added ØB

			# Convert all allele names in current marker to base pairs.
			basePairTmpLst[[marker]] <- offsetByColor[marker] + floor(alleleValue) * repeatUnitByColor[marker] + (alleleValue %% 1) * 10
		}
            bpmarkerByColorList[[color]] <- sapply(basePairTmpLst,function(x) x[1])
           
		# Add basepair to list.
		for(row in 1:length(allelesByColor)){

			# Avoid 'subscript out of bounds' error.
			if (length(basePairList)<color){
				basePairList[[color]] <- basePairTmpLst[[row]]
				phListByColor[[color]] <- peakHeightsByColor[[row]]
				allelesByColorList[[color]] <- allelesByColor[[row]]
                        if(isLab) refLabelByColorList[[color]] <- refLabelByColor[[row]]
			} else {
				basePairList[[color]] <- c(basePairList[[color]], basePairTmpLst[[row]])
				phListByColor[[color]] <- c(phListByColor[[color]], peakHeightsByColor[[row]])
				allelesByColorList[[color]] <- c(allelesByColorList[[color]], allelesByColor[[row]])
 				if(isLab) refLabelByColorList[[color]] <- c(refLabelByColorList[[color]], refLabelByColor[[row]])

			}
		}
	}


	# CREATE GRAPH

	# Set up the plot window according to the number of color channels.
	par(mfrow = c(noColors, 1))

	# Make x axis cross at y = 0 (i.e. remove the 4% margin at both ends of the y axis).
	par(yaxs = "i")

	# Reduce the spacing between the plots.
	# c(bottom, left, top, right) is the number of lines of margin to be specified on the four sides of the plot.
	# The default is c(5, 4, 4, 2) + 0.1
	par(mar = c(2.3, 4, 1.5, 1) + 0.1)

	# Define lower and upper bound for the x axis.
	xMin <- .Machine$integer.max 
	xMax <- 0
	for (row in 1:length(basePairList)){

#print("basePairList[[row]]")
#print(basePairList[[row]])

		xChannelMin <- sapply(na.omit(basePairList[[row]]),min)
		if (length(xChannelMin) > 0) {
			xChannelMin <- min(xChannelMin)
			xVal <- is.finite(xChannelMin)
		} else {
			xVal <- FALSE
		}
		if (xMin > xChannelMin && xVal == TRUE) {
			xMin <- xChannelMin
		}
		xChannelMax <- sapply(na.omit(basePairList[[row]]),max)
		if (length(xChannelMax) > 0) {
			xChannelMax <- max(xChannelMax)
			xVal <- is.finite(xChannelMax)
		} else {
			xVal <- FALSE
		}
		if (xMax < xChannelMax && xVal == TRUE) {
			xMax <- xChannelMax
		}
	}

	#Loop over all color channels.
	for (color in 1:noColors){

		# Get alleles and peak heights for current marker.
		bpVector <- basePairList[[color]]
		phList <- phListByColor[[color]]

		# Create blank plot with axes.
#		yMax <- sapply(na.omit(phList),max)
		yMax <- sapply(phList[!is.na(phList)],max)

		noData <- FALSE
		if (length(yMax) == 0) {
#print("length yMax == 0")
			yMax <- 1 # Minimal height.
			noData <- TRUE
		}
		yMax <- max(yMax)
#		yMax <- sapply(na.omit(yMax),max)

#		noData <- FALSE
		if (is.infinite(yMax) || is.na(yMax)) {
			yMax <- 1 # Minimal height.
			noData <- TRUE
		}
		plot(c(xMin, xMax), c(0, yMax), type="n", ylim = c(0, yMax * yMarginTop), ann = FALSE,axes=FALSE)
            axis(side = 2)
            nl <- 25 #number of ticks
            axis(side = 1,at=seq(xMin,xMax,l=nl),label=rep("",nl))
            if(isLab && color==1) legend("topright",legend=paste0("Label ",1:length(refnames)," = ",refnames),bty="n")
            abline(h=0)
            if(threshT>0) abline(h=threshT,col="gray",lwd=0.5)

		# Write text if no data.
		if (noData) {
			text(xMax / 1.4, yMax / 2, labels="No data", cex = 1.5)
		}

		# Create a title.
		if (color == 1) {
		title(main = sampleName, col.main = "red", font.main = 4)
		}	

		# Label the y axes.
		title(ylab = ylabel)

		# Label the x axis.
		# pos values of 1, 2, 3 and 4 indicate positions below, left, above and right of the coordinate.
#		mtext(paste(xlabel), side = 1, line = 2, adj = 0, cex = alleleNameTxtSize)
		mtext(paste(xlabel), side = 1, line = 0, adj = 0, cex = alleleNameTxtSize)

		# Write allele names under the alleles.
		# The additional par(xpd=TRUE) makes it possible to write text outside of the plot region.
		text(bpVector, 0, labels = allelesByColorList[[color]], cex = alleleNameTxtSize, pos = 1, xpd = TRUE) 
            if(isLab) text(bpVector,-yMax/20,labels=refLabelByColorList[[color]],cex = alleleNameTxtSize,pos=1, xpd = TRUE)

#           text(bpmarkerByColorList[[color]], yMax + yMax*0.02 ,expression(bold(paste0(markerByColorList[[color]]))),cex=alleleNameTxtSize,font=1)
            text(bpmarkerByColorList[[color]], yMax + yMax*0.04 ,markerByColorList[[color]],cex=1,font=2,xpd = TRUE)

		# Loop over all peaks.
		for (peak in 1:length(bpVector)){

			# Check if boxplot is to be drawn. Do not dra if only one value.
			if (drawBoxPlots && length(phList[[peak]]) > 1) {
			# Draw boxplots showing the distribution of peak heights.
			boxplot(phList[peak], add=TRUE, at=bpVector[peak],
				border=colors[color],pars=list(boxwex=peakHalfWidth*20, axes=FALSE))
			}


			if (drawPeaks) {
			# Create corners of peak triangle.
			xCords <- c(bpVector[peak] - peakHalfWidth, bpVector[peak], bpVector[peak] + peakHalfWidth)
			yCords <- c(0, lapply(phList[peak],mean), 0)
			# Plot peaks as filled polygons.
			polygon(xCords, yCords, col = tolower(as.character(colors[color])))
			}
		}
	}
 dev.new()
 op <- par(no.readonly = TRUE)
 dev.off()
 par(op)
} #end plot function

 #START FUNCTION
 alength <- sapply(Data,function(x) length(x$adata))
 hlength <- sapply(Data,function(x) length(x$hdata))
 locs <- names(Data)

 if( all(alength==0) ) {
  gmessage(message="There is no allele data in sample!",title="Error",icon="error")
 } else {

   #fix order prior:
   kit <- getKit(kitname)
   adata <- lapply(Data,function(x) x$adata)
   hdata <- lapply(Data,function(x) x$hdata)
   cdata <- list() #used to label conditioned references
   for(loc in locs) { #impute missing heights with 1
    lA <- length(adata[[loc]])
    if(length(hdata[[loc]])==0 && length(lA>0)) hdata[[loc]] <- rep(1,lA) #impute with height 1

    #reference data:
    if(!is.null(refcond)) {
     rdata <- lapply(refcond,function(x) unlist(x[[loc]]))
     unrdata <- unique(unlist(rdata))
     newA <- unrdata[!unrdata%in%adata[[loc]]] #get new alleles not seen in data
     nA <- length(newA)
     if(nA>0) {
      adata[[loc]] <- c(adata[[loc]], newA)
      hdata[[loc]] <- c(hdata[[loc]], rep(1e-6,nA)) #insert peak heights
     }
     cdata[[loc]] <- rep("",length(adata[[loc]]))
     for(rn in names(refcond)) {
       ri <- which(names(refcond)==rn) #reference index
       ind <- adata[[loc]]%in%rdata[[rn]]
       fix <- nchar(cdata[[loc]][ind])>0 
       cdata[[loc]][ind][fix] <- paste0(cdata[[loc]][ind][fix],"/")
       cdata[[loc]][ind] <- paste0(cdata[[loc]][ind],ri)
     }
    }
   } #end for each loci

   if (!is.na(kit[[1]][1])) { #the kit was found.
      kitlocs <- toupper(getKit(kitname,what="marker")) #make uppercase
      sname <- paste0(kitname," - ",sname)
      AMELname <- kitlocs[grep("AM",kitlocs,fixed=TRUE)] #get amel-name
      if(!is.null(AMELname)) {
       AMind <- grep("AM",toupper(locs),fixed=TRUE)
       locs[AMind] <- AMELname
       names(adata)[AMind] <- AMELname
       names(hdata)[AMind] <- AMELname
       if(!is.null(refcond)) names(cdata)[AMind] <- AMELname
      }
      #Insert missing loci:
      missloc <- kitlocs[!kitlocs%in%locs] #get missing loci
      for(loc in missloc) {
        adata[loc] = "" 
        hdata[loc] = 1e-6 #default is binary signal (below detection threshold)
        if(!is.null(refcond)) cdata[loc] = ""
      }
	adata <- adata[kitlocs] 
	hdata <- hdata[kitlocs] 
	if(!is.null(refcond)) cdata <- cdata[kitlocs] 
   } else {
    kitlocs <- locs
   }
   generateEPG(typingKit=kitname,alleleList=adata,peakHeightList=hdata, locusVector=kitlocs,sampleName=sname, drawBoxPlots=FALSE, drawPeaks=TRUE,refLabelList=cdata,refnames=names(refcond))
 }
}

