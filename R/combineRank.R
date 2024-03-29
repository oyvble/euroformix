#' @title combineRank
#' @author Oyvind Bleka
#' @description A search algorithm which creates a ranked table of the greatest combined sum from each vector-elements in a list.
#' @details Helpfunction used in earlier version of deconvolution to obtain a global ranking of full profiles (across all markers)

#' The procedure takes a list with vector elements of decreased sorted values and creates a ranked table of the most probable summed values (sum) from each vector-elements in the list.
#'
#' The search is stopped when SUM(exp(sum-loghdval))>alpha or when the size of the ranked table is equal maxsearch
#'
#' The search algorithm iterately compares a set new paths (by considering one step down in the ranked list for each list-elements). The information of earlier visited paths are stored with corresponing total sum. For the considered set of paths, the next iterative optimal path is the one with largest total sum.
#'
#' @param dlist A list with vector elements of decreased sorted values.
#' @param loghdval An optional value used to stop the search together with alpha.
#' @param alpha An optional value used to stop the search together with loghdval.
#' @param maxsearch The size of the ranked deconvolved profile table will not exceed this number (used to avoid endless search).
#' @return ret A list(rankG,pG) where rankG is the ranked genotypes with corresponding probabilities in pG.
#' @export
#' @keywords discrete optimalisation
combineRank <- function(dlist,loghdval=Inf,alpha=0.99,maxsearch=1000) {
 nD <- sapply(dlist,length)
 nL <- length(dlist)
 dmat <- matrix(-Inf,ncol=nL,nrow=max(nD))
 for(i in 1:nL) dmat[1:nD[i],i] <- dlist[[i]]

 #declare:
 trackmap <- numeric() #a matrix with index paths
 trackmapstring <- numeric() #a vector with stringed values of index paths
 trackval <- numeric() #a vector with sum to paths

 #init: 
 Rankmat <- rbind(rep(1,nL)) #recored rank matrix 
 sumdvec <- sum(dmat[1,]) #with corresponding sum score
 cc <- 2 #index counter

 while(cc<=maxsearch) {
  step <- rep(1,nL)
  deadend <- nD<(Rankmat[cc-1,] + step) #check if we have reached the dead end for some nodes (they are excluded as possible paths)
  if(all(deadend)) break #stop if no more paths to go!

  nextpaths <-  matrix(rep(Rankmat[cc-1,],nL),nrow=nL,byrow=TRUE) + diag(step) #see forward paths
  nextpaths <- nextpaths[!deadend,] #remove dead ends
  if(sum(!deadend)==1)  nextpaths <- t(nextpaths) #make matrix again

  #Store new possible tracks with values to map: (but check if nextpaths is not already in trackmap)
  for(k in 1:nrow(nextpaths)) {
   check <- paste0(nextpaths[k,],collapse="-")
   if(check%in%trackmapstring) next; #the path was already saved in the trackmap. We don't save it again!
   trackmapstring <- c(trackmapstring,check) #add string-key
   trackval <- c(trackval,sum(dmat[cbind(nextpaths[k,],1:nL)])) #add sum of new paths 
   trackmap <- rbind(trackmap,nextpaths[k,]) #add path to map
  } 
  #Select best path from map:
  pathsel <- which.max(trackval)[1] #could be multiple paths with same values!
  Rankmat <- rbind(Rankmat, trackmap[pathsel,]) #select path
  sumdvec <- c(sumdvec,trackval[pathsel]) #store corresponding value

  if(sum(exp(sumdvec-loghdval))>=alpha) break #searching criterion met!

  #remove selected path from map:
  trackmap <- trackmap[-pathsel,] #remove selected path from map
  trackval <- trackval[-pathsel] #remove corresponding value path from map
  trackmapstring  <- trackmapstring[-pathsel]
  cc <- cc +1  #go to next
 }
 return(list(rankG=Rankmat,pG=exp(sumdvec-loghdval)))
}
