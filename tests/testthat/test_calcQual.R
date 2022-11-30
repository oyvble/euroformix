#Test that calcQual calculates as expected (comparison with R-script)
#rm(list=ls())

#nUnknowns=2;nAlleles=2; nConds=0; nKnowns=0; nReps=2; pD=0.1; pC=0.05; fst=0.01
#nUnknowns=0;nAlleles=2; nConds=1; nKnowns=0; nReps=1; pD=0.1; pC=0.05; fst=0.01
compare = function(nUnknowns=4, nAlleles = 6, nConds = 0, nKnowns = 0, nReps=1, pD=0.1, pC=0.05, fst=0.01 ) {
  set.seed(1)
  evid = sample(2:30,nAlleles,replace = FALSE) #1:nA
  freq = setNames(rep(1/(nAlleles+1),nAlleles+1), c(evid,"99"))
  if(nAlleles==0) evid=0
  
  evids = evid 
  #TESTING SPECIAL CASE OF REPLICATES
  if(nReps>1) {
    for(rr in 2:nReps) {
      evids = c(evids,0,evid)
    }
  } else if(nReps<1) { #test special case:
    #First evid as normal, last is missing
    evids = c(evids,0,0)
    if(nReps==-1) {
      evids = c(evids,0,evid) #add another replicate
    }
  }
  
  conds = NULL # evids[3:4]
  knowns = NULL #evids[1:2]
  if(nConds>0)  conds = sample(names(freq),2*nConds,replace = TRUE)
  if(nKnowns>0)  knowns = sample(names(freq),2*nKnowns,replace = TRUE)

  #FOrensim first:
  nContr = length(conds)/2 + nUnknowns
  pDvec = rep(pD,nContr)
  #lik0 = forensim::likEvid(evids,conds,knowns,nUnknowns,fst, pDvec, pDvec^2,pC,freq) 
  lik0 = likQualR(evids,conds,knowns,nUnknowns,fst, pD, pC, freq) 
  lik1 = calcQual(evids,conds,knowns,nUnknowns,fst, pD, pC, freq) 
  #lik2 = euroformix::calcQual(evids,conds, knowns, nUnknowns, fst, pDvec, pC, freq) 
  return(c(lik0,lik1))#,abs(lik0-lik2)))
}  
#compare(nAlleles=2, nReps=2)

#TRAVERSIN FOLLOWING VECTORS (ALL COMBINATINOS)
pDV = c(0,0.1)
pCV = c(0,0.05)#,1)
fstV = c(0,0.01)#,1)
nCondsV = 0:2
nKnownsV = 0:1
nUnknownsV = 0:2
nAllelesV = c(0,1,4)
nRepsV = -1:2

test_that("C-implementation of qualitative model is correct", {
  
for(a in seq_along(nUnknownsV)) {
  for(b in seq_along(nAllelesV)) {
    for(c in seq_along(nCondsV)) {
      for(d in seq_along(nKnownsV)) {
        for(e in seq_along(nRepsV)) {
          for(f in seq_along(pDV)) {
            for(g in seq_along(pCV)) {
              for(h in seq_along(fstV)) {
                nUnknowns=nUnknownsV[a]
                nAlleles=nAllelesV[b]
                nConds=nCondsV[c]
                nKnowns=nKnownsV[d]
                nReps=nRepsV[e]
                pD=pDV[f]
                pC=pCV[g]
                fst=fstV[h]
                vals = compare(nUnknowns, nAlleles, nConds, nKnowns, nReps, pD, pC, fst)
                expect(vals[1],vals[2],10)
              }
            }
          }
        }
      }
    }
  }
}
  
})
#View(tab)

