#Testing that the numerical inference of gamma distribution is working
#rm(list=ls());library(euroformix);library(testthat)


th = c(1000,0.8,0.6) #true parameters
n = 30 #number of samples
x = seq(10,300,l=n)
#set.seed(1);y=rgamma(n,shape=(2/th[2]^2)*th[3]^((x-125)/100),scale=th[1]*th[2]^2)
y = c(2434.10045389206,5255.45799208187,4930.54859484562,3346.79665881895,4516.40953253601,3482.67267894867,3084.75989033535,1849.68456609539,1253.68266225878,1319.74306160061,1528.96067970406,1684.7590275416,845.986838779711,2453.67407627856,2066.84798535382,2342.70302678474,2072.08528478055,1256.09166426717,67.6821918881922,1495.59122820608,2022.12297997056,459.826219293013,533.707145276646,632.617057156043,667.27247237923,321.461122468519,1241.59025314987,352.539656179777,1381.81288737032,1022.15471305975)

test_that("fitted gamma distr with degrad :", {
  thhat1 = fitgammamodel(y,x,DEG=TRUE,niter=1)#,offset = 125, scale = 100,plott=T)
  expect_equal(thhat1,c(968.9315521 ,  0.6820142  , 0.5135972))
})

test_that("fitted gamma distr without degrad :", {
  thhat2 = fitgammamodel(y,x,DEG=FALSE,niter=1)#,offset = 125, scale = 100,plott=T)
  expect_equal(thhat2,c(932.052467, 1.080866))
})

