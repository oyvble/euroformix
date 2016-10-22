#' @title getKit
#' @author Oskar Hansson
#' @description Function to get kit information. Tided up by Øyvind Bleka
#' @details Returns kit information
#' @param kit shortname of kit: {"ESX17","ESI17","ESI17Fast","ESX17Fast","Y23","Identifiler","NGM","ESSPlex","ESSplexSE","NGMSElect","SGMPlus","ESX16", "Fusion","GlobalFiler"}
#' @param what options: {"Index","Panel","Short.Name","Full.Name","Marker","Allele","Size","Virtual","Color","Repeat","Range","Offset","Gender"}
#' @return res A data frame with kit information
#' @export

getKit <- function(kit=NULL, what=NA) {  
  .separator <- .Platform$file.sep # Platform dependent path separator. 
 #  packagePath <- "C:/Users/oebl/Dropbox/Forensic/MixtureProj/myDev/quantLR/euroformix0"
  packagePath <- path.package("euroformix", quiet = FALSE) # Get package path.
  subFolder <- "extdata"
  fileName <- "kit.txt"
  filePath <- paste(packagePath, subFolder, fileName, sep=.separator)
  .kitInfo <- read.delim(file=filePath, header = TRUE, sep = "\t", quote = "\"",dec = ".", fill = TRUE, stringsAsFactors=FALSE)
 
  # Available kits. Must match else if construct.
  kits<-unique(.kitInfo$Short.Name)
	if (is.null(kit)) {	# Check if NULL
		res<-kits
	} else {	# String provided.
		# Check if number or string.
		if (is.numeric(kit)) {
			index<-kit # Set index to number.
		} else {
			index<-match(toupper(kit),toupper(kits)) # Find matching kit index (case insensitive)
		}
		if (any(is.na(index))) { 		# No matching kit.
			return(NA)
		# Assign matching kit information.
		} else {
		  currentKit <- .kitInfo[.kitInfo$Short.Name==kits[index], ]
              res <- data.frame(Panel = currentKit$Panel,
                        Short.Name = currentKit$Short.Name,
                        Full.Name = currentKit$Full.Name,
                        Marker = currentKit$Marker,
                        Allele = currentKit$Allele,
                        Size = currentKit$Size,
                        Size.Min = currentKit$Size.Min,
                        Size.Max = currentKit$Size.Max,
                        Virtual = currentKit$Virtual,
                        Color = currentKit$Color,
                        Repeat = currentKit$Repeat,
                        Marker.Min = currentKit$Marker.Min,
                        Marker.Max = currentKit$Marker.Max,
                        Offset = currentKit$Offset,
                        Gender.Marker = currentKit$Gender.Marker,
                        stringsAsFactors = FALSE)
		  res$Marker <- factor(res$Marker, levels=unique(res$Marker))# Create useful factors. 
		} 
	}
      
  # Kit is required.
  if (!is.null(kit)) {

    if(is.na(what)){  # Return all kit information.
      return(res)
    } else if (toupper(what) == "INDEX"){ # Return kit index.
      return(index)
    } else if (toupper(what) == "PANEL"){  # Return panel name.
      return(unique(res$Panel))
    } else if (toupper(what) == "SHORT.NAME"){ # Return short name.
      return(unique(res$Short.Name))
    } else if (toupper(what) == "FULL.NAME"){ # Return full name.
      return(unique(res$Full.Name))
    } else if (toupper(what) == "MARKER"){  # Return all markers.
      return(as.vector(unique(res$Marker)))
    } else if (toupper(what) == "ALLELE"){ # Return all alleles and markers. 
      res <- data.frame(Marker=res$Marker, Allele=res$Allele)
      return(res)
    } else if (toupper(what) == "SIZE"){
      # Returns all alleles and their indicated normal size in base pair.
      # Their normal size range is idicated in min and max columns.
      # Grouped by marker.
      res <- data.frame(Marker=res$Marker,Allele=res$Allele,Size=res$Size, Size.Min=res$Size.Min, Size.Max=res$Size.Max,stringsAsFactors=FALSE)
      return(res)
    } else if (toupper(what) == "VIRTUAL"){
      # Returns all alleles (bins) with a flag if it is virtual
      # 1 for virtual or 0 it it is a physical ladder fragment.
      # Grouped per marker.
      res <- data.frame(Marker=as.character(res$Marker), Allele=res$Allele,Virtual=res$Virtual, stringsAsFactors=FALSE)
      return(res)
    } else if (toupper(what) == "COLOR"){ # Return markers and their color as strings.
      marker <- getKit(kit, what="Marker")
      color <- NA

      for(m in seq(along=marker)){
        color[m] <- unique(res$Color[res$Marker == marker[m]])
      }
      res <- data.frame(Marker=marker,Color=color,stringsAsFactors=FALSE)
      return(res)
      
    } else if (toupper(what) == "REPEAT"){ # Return markers and their repeat unit length in base pair.
      marker <- getKit(kit, what="Marker")
      offset <- NA
      repeatUnit <- NA
      
      for(m in seq(along=marker)){
        offset[m] <- unique(res$Offset[res$Marker == marker[m]])
        repeatUnit[m] <- unique(res$Repeat[res$Marker == marker[m]])
      }
      res <- data.frame(Marker=marker, Offset=offset, Repeat=repeatUnit, stringsAsFactors=FALSE)
      return(res)
      
    } else if (toupper(what) == "RANGE"){   # Return markers and their range (min and max) in base pair.
      marker <- getKit(kit, what="Marker")
      markerMin <- NA
      markerMax <- NA
      color <- NA
    
      for(m in seq(along=marker)){
        markerMin[m] <- unique(res$Marker.Min[res$Marker == marker[m]])
        markerMax[m] <- unique(res$Marker.Max[res$Marker == marker[m]])
        color[m] <- unique(res$Color[res$Marker == marker[m]])
      }
      res <- data.frame(Marker=marker, Color=color,  Marker.Min=markerMin,  Marker.Max=markerMax,stringsAsFactors=FALSE)
      res$Color <- factor(res$Color, levels=unique(res$Color)) # Create useful factors.
      return(res)
      
    } else if (toupper(what) == "OFFSET"){  # Return markers and their estimated offset in base pair.
      marker <- getKit(kit, what="Marker")
      offset <- NA
      repeatUnit <- NA
  
      for(m in seq(along=marker)){
        offset[m] <- unique(res$Offset[res$Marker == marker[m]])
        repeatUnit[m] <- unique(res$Repeat[res$Marker == marker[m]])
      }
      res <- data.frame(Marker=marker, Offset=offset, Repeat=repeatUnit,stringsAsFactors=FALSE)
      return(res)
      
    } else if (toupper(what) == "GENDER"){  # Return gender marker as string. 
      genderMarker <- as.character(unique(res$Marker[res$Gender.Marker == TRUE]))
      if(length(genderMarker) > 1){
        warning(paste("More than one gender marker returned for kit", kit))
      }
      return(genderMarker)
    } else {
      warning(paste(what, "not supported! \nwhat = {", options,"}"))
      return(NA)
	}
 } else { # If kit is NULL return available kits.
   return(res)  
 }
} #end function

