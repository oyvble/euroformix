#' @title setVecRightOrder
#' @author Oyvind Bleka
#' @description Sorting vector with respect to names of an input (case insensitive)
#' @details Helpfunction to set vector to same locus order as output from prepareC (order based on evidence profiles)
#' Used in functions contLikMLE, contLikINT, calcRMPfst
#' @param vec A vector with values to sort by letting names(vec) have same order as orderedNames
#' @param orderedNames A vector with names to sort vec vector with
#' @return Vector with same order as specified in orderedNames
#' @export

setVecRightOrder = function(vec,orderedNames) { 
  if(is.null(names(vec))) stop("The setting vector must must contain locus names (any order)!") #stop if locus names not given
  vec = vec[toupper(names(vec))%in%toupper(orderedNames)] #keep only names specified in orderedNames
  if(length(vec)!=length(orderedNames)) stop("The setting vector was not same length as the number of markers (mismatch)!")
  ord = match(toupper(orderedNames),toupper(names(vec)))
  if(any(is.na(ord))) stop("The setting vector was missing locus name. Please check!")
  return( vec[ord]) #get vector with correct order
}