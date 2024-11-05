#' @title prepareCobj
#' @author Oyvind Bleka
#' @description Wrapper function to prepare the C++ module 
#' @details Takes values from prepareC functions and prepares datastructure in C++
#' @param c returned list object from prepareC function
#' @param verbose Whether to printing time used
#' @return Pointer to C++ module
#' @export 

prepareCobj = function(c,verbose=FALSE) {
  mod = Rcpp::Module( "mod",PACKAGE="euroformix" ) #load module
  obj = methods::new(mod$ExposedClass) #create object of class
  obj$filldata(c$nStutterModels,c$nMarkers,c$nRepMarkers,c$nAlleles,c$startIndMarker_nAlleles,c$startIndMarker_nAllelesReps,c$peaks,c$freqs,c$dropinWeight, c$nTyped, c$maTyped, c$basepair,
               c$BWfrom, c$FWfrom, c$BWto, c$FWto, c$nPotStutters, c$startIndMarker_nAllelesTot, c$QalleleIndex, c$dropinProb, c$fst, c$AT, c$lambda, c$NOK, c$knownGind, c$relGind, c$ibd, c$maxThreads) 
  prepTime=system.time({
    obj$prepare(as.integer(c$nC))
  })[3]
  if(verbose) print(paste0("Prep. done and took ",round(prepTime),"s. Start pre-search..."))
  return(obj)
} #end function

