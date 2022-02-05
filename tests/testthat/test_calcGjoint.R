#rm(list=ls());library(euroformix);library(testthat)

fst=0.1 #use high value
freq = c(0.1,0.2,0.7 )
names(freq) = seq_along(freq)
refK =  c(1,2,3,3) #typed
refR = c(1,2) #ref alleles
dec=4 #round off error

#PART 1: #one unknown
nU = 1

#paste0(round(x,dec),collapse=",")
test_that("Unrelated(nU=1)", {
  ibd = c(1,0,0) 
  Glist = calcGjoint(freq, nU, fst, refK, refR, ibd) #include all typed refs here
  x = round(Glist$Gprob,dec) 
  expect_equal(  as.numeric(x), c(0.0303,0.0585,0.1733,0.0585,0.2554,0.4241))
})
test_that("Sibling(nU=1)", {
  ibd = c(0.25,0.5,0.25) 
  Glist = calcGjoint(freq, nU, fst, refK, refR, ibd) #include all typed refs here
  x = round(Glist$Gprob,dec) 
  expect_equal(  as.numeric(x), c(0.0441,0.355,0.2029,0.0685,0.2235,0.106))
})
test_that("Child(nU=1)", {
  ibd = c(0,1,0) 
  Glist = calcGjoint(freq, nU, fst, refK, refR, ibd) #include all typed refs here
  x = round(Glist$Gprob,dec) 
  expect_equal(  as.numeric(x), c(0.0731,0.1808,0.3192,0.1077,0.3192,0))
})

#PART 2: #two unknowns
nU = 2 
#paste0(round(c(x),dec),collapse=",")
test_that("Unrelated(nU=2)", {
  ibd = c(1,0,0) 
  Glist = calcGjoint(freq, nU, fst, refK, refR, ibd) #include all typed refs here
  x = round(Glist$Gprob,dec) 
  val = c(0.0024,0.0028,0.0082,0.0013,0.0059,0.0097,0.0028,0.0054,0.0117,0.0044,0.0154,0.0188,0.0082,0.0117,0.0389,0.0077,0.0376,0.0692,0.0013,0.0044,0.0077,0.0068,0.0194,0.0188,0.0059,0.0154,0.0376,0.0194,0.0752,0.1019,0.0097,0.0188,0.0692,0.0188,0.1019,0.2057)
  expect_equal(  as.numeric(x), val)
})
test_that("Sibling(nU=2)", {
  ibd = c(0.25,0.5,0.25) 
  Glist = calcGjoint(freq, nU, fst, refK, refR, ibd) #include all typed refs here
  x = round(Glist$Gprob,dec) 
  val = c(0.0026,0.0035,0.0104,0.0022,0.0096,0.0159,0.0116,0.0225,0.0627,0.0223,0.0919,0.144,0.0062,0.011,0.0366,0.01,0.049,0.0901,0.0017,0.0048,0.01,0.0064,0.021,0.0245,0.0057,0.0119,0.0363,0.0129,0.0584,0.0983,0.0024,0.0047,0.0173,0.0047,0.0255,0.0514)
  expect_equal(  as.numeric(x), val)
})
test_that("Child(nU=2)", {
  ibd = c(0,1,0) 
  Glist = calcGjoint(freq, nU, fst, refK, refR, ibd) #include all typed refs here
  x = round(Glist$Gprob,dec) 
  val = c(0.0039,0.0057,0.0168,0.0037,0.0162,0.0269,0.0068,0.0131,0.0329,0.0131,0.0485,0.0664,0.0084,0.0162,0.0537,0.0162,0.0792,0.1456,0.0028,0.0074,0.0162,0.0094,0.0323,0.0396,0.0084,0.0162,0.0537,0.0162,0.0792,0.1456,0,0,0,0,0,0)
  expect_equal(  as.numeric(x), val)
})


