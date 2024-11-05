#' @title plotTopEPG
#' @author Oyvind Bleka
#' @description Replaced by plotTopEPG2
#' @param MLEobj An object returned from calcMLE/contLikMLE
#' @param DCobj An object returned from devonvolve: Must be run with same object as MLEobj
#' @param kitname Name of kit:  Obtained from getKit()
#' @param threshT The detection threshold can be shown in gray in the plot.
#' @export

plotTopEPG <- function(MLEobj,DCobj=NULL,kitname=NULL,threshT=0) {
  plotTopEPG2(MLEobj,DCobj,kitname) #run this function instead
}

