#
#' @title getMarkerVal 
#' @author Oyvind Bleka
#' @description Extract value with respect to names of an input (case insensitive)
#' @details Helpfunction to extract a specific value from a vector wrt a specific locus 
#' Used in functions genDataset, prepareData
#' @param vec A vector with values to sort by letting names(vec) have same order as orderedNames
#' @param marker The marker to extract value of
#' @return the value for selected marker
#' @export

getMarkerVal = function(vec,marker="") { #helpfunction to obtain marker based value
  errmsg = "The argument must contain valid locus name in names!"
  if(length(vec)==1) { #if exactly one element
    return(vec) #value is fine
  } else if(length(vec)==0) {
    stop("No value found in argument!")
  } else { #if more than 1 element
    
    if(is.null(names(vec))) {
      stop(errmsg)
      
    } else { #if names was included
      x0 = vec[toupper(names(vec))==toupper(marker)] #obtain value for corresponding marker
      if(length(x0)==1) {
        return(x0) #return value
      } else {
        stop(errmsg)
      }
    }
  }
} #end funciton
