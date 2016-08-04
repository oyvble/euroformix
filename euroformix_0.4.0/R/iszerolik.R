#' @title iszerolik 
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description iszerolik Determine if the likelihood becomes zero
#' @details The function is determining whether likelihood vil be 0 (TRUE or FALSE) assuming zero drop-in 
#' @param evid A vector with allele names
#' @param ref A vector with conditioned alleles in given hypothesis
#' @param nU Number of unknown individuals in given hypothesis
#' @return TRUE/FALSE A boolean whether the likelihood vil be zero.
#' @export 

iszerolik <- function(evid,ref,nU,xi=0) {
   Ei2 <- evid[!evid%in%ref] #set of unknown alleles (not explained by ref0) 
   if(length(Ei2)<=(2*nU)) return(FALSE) #unknown contributors explains E
   if(!is.null(xi) && xi==0) return(TRUE) #too many unexplained alleles
   #case of assuming stutters:
   Ei3 <- Ei2[!as.character(Ei2)%in%as.character(as.numeric(ref)-1)] #set of unknown alleles (after explained by being stutter from ref0)
   avec <- sort(as.numeric(Ei3),decreasing=TRUE) #sort allele names in decreasing order
   c <- 0 #counter of number of alleles required to explain with stutter:
   while(1) {
    if(length(avec)==0) break #stop loop when done
    a <- avec[1]
    avec <- avec[!as.character(avec)%in%as.character(c(a,a-1))] #sequential remove alleles explainable as stutters
    c <- c+1
    if(c>(2*nU)) return(TRUE)  #too many unexplained alleles
   }
   return(FALSE) #the alleles can be explained
}
