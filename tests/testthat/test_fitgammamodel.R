#Testing that the numerical inference of gamma distribution is working
#NEW VERSION FOR v4.1.0
#rm(list=ls());library(euroformix);library(testthat)

#PART 1: beta < 1
n = 30 #number of samples
x = seq(10,300,l=n)
dec = 3 #decimals numbers

set.seed(1)
th = c(1000,0.8,0.6) #true parameters
y=rgamma(n,shape=(2/th[2]^2)*th[3]^((x-125)/100),scale=th[1]*th[2]^2)
#y = c(2434.10045389206,5255.45799208187,4930.54859484562,3346.79665881895,4516.40953253601,3482.67267894867,3084.75989033535,1849.68456609539,1253.68266225878,1319.74306160061,1528.96067970406,1684.7590275416,845.986838779711,2453.67407627856,2066.84798535382,2342.70302678474,2072.08528478055,1256.09166426717,67.6821918881922,1495.59122820608,2022.12297997056,459.826219293013,533.707145276646,632.617057156043,667.27247237923,321.461122468519,1241.59025314987,352.539656179777,1381.81288737032,1022.15471305975)

set.seed(1)
test_that("fitted gamma distr with degrad:", {
  thhat1 = fitgammamodel(y,x)#,offset = 125, scale = 100,plott=T)
  expect_equal(round(thhat1,dec),c(968.935 ,  0.682  , 0.514))
})

set.seed(1)
test_that("fitted gamma distr without degrad :", {
  thhat2 = fitgammamodel(y)#,offset = 125, scale = 100,plott=T)
  expect_equal(round(thhat2,dec),c(932.056, 1.081))
})


#PART 2: beta > 1
th = c(1000,0.8,1.2) #true parameters
set.seed(1)
y=rgamma(n,shape=(2/th[2]^2)*th[3]^((x-125)/100),scale=th[1]*th[2]^2)
#plot(x,y)
set.seed(1)
test_that("fitted gamma distr with degrad (noDEG) :", {
  thhat1 = fitgammamodel(y,x,restrictDeg = TRUE)#,offset = 125, scale = 100,plott=T)
  expect_equal(round(thhat1,dec),c(1007.097  ,  0.707  , 0.999))
})

set.seed(1)
test_that("fitted gamma distr with degrad (noDEG):", {
  thhat2 = fitgammamodel(y,x,restrictDeg = FALSE)#,offset = 125, scale = 100,plott=T)
  expect_equal(round(thhat2,dec),c(1007.097,    0.707, 1.059))
})
