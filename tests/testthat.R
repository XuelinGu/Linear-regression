library(testthat)
library(stats)
library(Linearregression)
y=c(1,2,3,4,5)
x1=c(1,2,3,4,5)
x2=c(3,4,2,5,2)

test_that("Linearregression", {
  expect_equal(lr(y~x1+x2)[[1]][1,1],lm(y~x1+x2)[[1]][[1]])
})
