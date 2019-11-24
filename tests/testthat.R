library(testthat)
library(stats)
library(Linearregression)

####testing small samples
y = c(23,24,26,37,38,25,36,40)
x1 = c(1,2,3,4,5,6,7,8)
x2 = c(23,32,34,20,24,56,34,24)
x3 = c("M","F","F","U","M","F","U","M")
x4 = mtcars$hp
test_that("Linearregression", {
  ##compare the estimated coeffecients
  expect_equal(as.numeric(lr(y~x1+x2)[[1]][,1]),as.numeric(lm(y~x1+x2)[[1]]))
  expect_equal(as.numeric(lr(y~x1+x2+x3,reference='F')[[1]][,1]),as.numeric(lm(y~x1+x2+x3)[[1]]))
  expect_equal(as.numeric(lr(mtcars$mpg~mtcars$cyl+mtcars$disp)[[1]][,1]),as.numeric(lm(mtcars$mpg~mtcars$cyl+mtcars$disp)[[1]]))
  expect_equal(as.numeric(lr(mpg~cyl+disp+x4,data = mtcars)[[1]][,1]),as.numeric(lm(mpg~cyl+disp+x4,data = mtcars)[[1]]))
})



####testing big samples
set.seed(4)
y = rnorm(1e+3,1,10)
x1 = rnorm(1e+3,1,10)
x2 = rnorm(1e+3,1,10)
test_that("Linearregression", {
  ##compare the estimated coeffecients
  expect_equal(as.numeric(lr(y~x1+x2)[[1]][,1]),as.numeric(lm(y~x1+x2)[[1]]))
  })

