#' @title plotTopLUS
#' @author Oyvind Bleka
#' @description A dataplotter for MPS data with LUS format
#' @details Utilizing existing function plotTopMPS2 to show.
#' @param MLEobj An object returned from calcMLE/contLikMLE
#' @param DCobj An object returned from deconvolve: Must be run with same object as MLEobj
#' @param threshT The detection threshold can be shown in gray in the plot.
#' @param LUSsymbol The separator for each variant. First is assumed to be CE.
#' @export

plotTopLUS = function(MLEobj,DCobj=NULL,threshT=0,LUSsymbol="_") {
  plotTopMPS2(MLEobj,DCobj,grpsymbol=LUSsymbol) #simple call
} #end function

