#library(euroformix)
#source("C:/Users/oyvbl/Dropbox/Forensic/MixtureProj/myDev/euroformix/R/helpintegrals.R", echo=TRUE)

calcMyInt = function(nC,nKnown=0) {
  
#  nC = 3;nKnown=0
  #Init an own integration function
  myInt = function(fun,L,U,...) {
    int = cubature::adaptIntegrate(fun,L,U, tol=1,...)[[1]]
    #int = integrate(Vectorize(fun),L,U, rel.tol=1,...)[[1]]
    return(int)
  }

  #
  #DEFINING THE LIKELIHOOD EXPRESSION: dir(alpha=1). 
  alphaVec=rep(1,nC) #any choice is ok
  logConst = lgamma(sum(alphaVec)) - sum(lgamma(alphaVec))
  const = exp(logConst)
  nEvals = 0
  LikFun = function(par) { #helpfunction to obtain (Scaled) likelihood value
    nEvals <<- nEvals + 1
    par2 = c(par,1-sum(par))
    #val = exp(sum((alphaVec-1)*log(par2)) + logConst)
    val = const
    return(val)
  }
  
  mxLims = NULL
  if(nC>1) { #only relevant for having at least 2 contributors
    mxLims = cbind(0,1/(2:nC)) #obtain Mx limits in case of unknowns
    if(nKnown>0) { #if conditionals
      mxLimsKnown = t(replicate(nKnown,c(0,1))) #get limits for knowns
      mxLims = rbind(mxLimsKnown,mxLims)[1:(nC-1),,drop=FALSE]
    }
  }

  #Indicate calculation situations
  calcIntMx = NULL
  if(nC==1) calcIntMx = euroformix:::.calcInt1p
  if(nC==2) {
    if(nKnown==0) calcIntMx = euroformix:::.calcInt2p0K
    if(nKnown>0) calcIntMx = euroformix:::.calcInt2p1K
  }
  if(nC==3) {
    if(nKnown==0) calcIntMx = euroformix:::.calcInt3p0K
    if(nKnown==1) calcIntMx = euroformix:::.calcInt3p1K
    if(nKnown>1) calcIntMx = euroformix:::.calcInt3p2K
  }
  if(nC==4) {
    if(nKnown==0) calcIntMx = euroformix:::.calcInt4p0K
    if(nKnown==1) calcIntMx = euroformix:::.calcInt4p1K
    if(nKnown==2) calcIntMx = euroformix:::.calcInt4p2K
    if(nKnown>2) calcIntMx = euroformix:::.calcInt4p3K
  }
  if(nC==5) {
    if(nKnown==0) calcIntMx = euroformix:::.calcInt5p0K
    if(nKnown==1) calcIntMx = euroformix:::.calcInt5p1K
    if(nKnown==2) calcIntMx = euroformix:::.calcInt5p2K
    if(nKnown==3) calcIntMx = euroformix:::.calcInt5p3K
    if(nKnown>3) calcIntMx = euroformix:::.calcInt5p4K
  }

  int = calcIntMx(LikFun,myInt,mxLims,thetaOther = NULL)
  return(list(int=int,nEvals=nEvals))  
}

test_that("Integrals are defined as expected", {
  for(NOC in 1:5) {
    for(nCond in 0:NOC) {
      #    NOC=5;nCond=3
      #print(paste0("NOC=",NOC," - nCond=",nCond))
      ret = calcMyInt(NOC,nCond)
      expect_equal(ret$int,1)
    }
  }
})

