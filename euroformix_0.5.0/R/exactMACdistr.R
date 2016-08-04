#' @title exactMACdistr
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description Calculation of false positive probability of Matching allele counter.
#' @details This function will calculate the false positive probability for a given range of allele matches.
#' Thanks 1: Thore Egeland: Explicit formulas.
#' Thanks 2: Aaron (stat.umn.edu/~arendahl): A most useful quick permutation-implementation
#'
#' @param si A vector with each element as the sum of allele frequencies for observed alleles.
#' @param lowk Lower range of number of allele matches.
#' @return ret A vector with false positive probabilities for each given number of observed alleles.
#' @export


exactMACdistr <- function(si,lowk=0) {
 #si - A vector: With elements equal the sum of the observed frequencies in Evidence

 #function for calculing exact (n-k) false positive MAC
 exactMAC = function(si,k) {
  #si - sum of observed frequencies used in calculation
  #k - number of missmatching alles
  #Recursive method returning sequences:
   #unique permutations: found at:
   #"http://stackoverflow.com/questions/5671149/permute-all-unique-enumerations-of-a-vector-in-r"
   uniqueperm2 <- function(d) {
    dat <- factor(d)
    N <- length(dat)
    n <- tabulate(dat)
    ng <- length(n)
    if(ng==1) return(d)
    a <- N-c(0,cumsum(n))[-(ng+1)]
    foo <- lapply(1:ng, function(i) matrix(combn(a[i],n[i]),nrow=n[i]))
    out <- matrix(NA, nrow=N, ncol=prod(sapply(foo, ncol)))
    xxx <- c(0,cumsum(sapply(foo, nrow)))
    xxx <- cbind(xxx[-length(xxx)]+1, xxx[-1])
    miss <- matrix(1:N,ncol=1)
    for(i in seq_len(length(foo)-1)) {
      l1 <- foo[[i]]
      nn <- ncol(miss)
      miss <- matrix(rep(miss, ncol(l1)), nrow=nrow(miss))
      k <- (rep(0:(ncol(miss)-1), each=nrow(l1)))*nrow(miss) + l1[,rep(1:ncol(l1), each=nn)]
      out[xxx[i,1]:xxx[i,2],] <- matrix(miss[k], ncol=ncol(miss))
      miss <- matrix(miss[-k], ncol=ncol(miss))
    }
    k <- length(foo)
    out[xxx[k,1]:xxx[k,2],] <- miss
    out <- out[rank(as.numeric(dat), ties="first"),]
    foo <- cbind(as.vector(out), as.vector(col(out)))
    out[foo] <- d
    t(out)
  }
  recSequence = function(x,S) {
   #x working sequence with (1,2)
   #S is total sum
   #X is stored sequences
   if(S==0) {
    xList <<- list(numeric())
    return(NULL) #no sequences satisfied
   }
   if(sum(x)<S) { #still more to go
    for(i in 1:2) { #adding possibilities
     xList <- recSequence( c(x,i),S)
    } 
   } else if(sum(x)==S) { #accept sequence
    xList[[length(xList)+1]] <<- x
    return
   } else { #don't accept sequence
    return
   }
  } 
  n <- length(si) #number of loci
  pZ0 <- (1-si)^2
  pZ1 <- 2*si*(1-si)
  pZ2 <- si^2 
  xList <<- list()
  recSequence(x=numeric(),k)

  #need only unique sequences and <=n:
  tmp <- numeric()
  yList <- list()
  for(i in 1:length(xList)) {
   L <- length(xList[[i]])
   if(L <= n && !any(tmp==L)) {
    tmp <- c(tmp,L)
    yList[[length(tmp)]] <- xList[[i]]
   }
  }
  SS <- numeric()
  for(i in 1:length(yList)) {
   x <- yList[[i]]
   Lx <- length(x)
   Ly <- n-Lx
   y <- rep(0,Ly)
   Z <- uniqueperm2(c(x,y)) #add together to all possible combinations
   if(is.null(dim(Z))) Z <- t(Z)
   W <- matrix(NA,nrow=nrow(Z),ncol=ncol(Z)) 
   for(j in 1:n) {
    Zi <- Z[,j]
    W[Zi==2,j] <- pZ0[j] #no alleles
    W[Zi==1,j] <- pZ1[j] #one allele
    W[Zi==0,j] <- pZ2[j] #both alles
   }
   w <- exp(rowSums(log(W))) #product over columns
   SS <- c(SS,sum(w))
  }
  return(sum(SS))
 }
 
 #function starts:
 n <- length(si)
 if(lowk>(2*n) | lowk<0) print("Specify lowk to be between 0 and 2*#markers")
 kvec <- lowk:(2*n) #range of what we want
 prob_k <- rep(NA,length(kvec))
 for(k in 1:length(kvec)) {
  prob_k[k] <- exactMAC(si,(2*n)-kvec[k])
 } 
 names(prob_k) <- kvec #number of matches
 return(prob_k)
}



