library(testthat)
library(stats)
library(Linearregression)

##lm() in stats package is commonly used for linear regression fitting so we could use it to test the accuracy of lr().
##Since a list of results were produced both in lr() and lm() , we only test accuracy of the estimated coeffecients.

####1 small samples
y = c(23, 24, 26, 37, 38, 25, 36, 40)
x1 = c(1, 2, 3, 4, 5, 6, 7, 8)
x2 = c(23, 32, 34, 20, 24, 56, 34, 24)
x3 = c("M", "F", "F", "U", "M", "F", "U", "M")
x4 = mtcars$hp

##1.1 cope with different compatible formula format
test_that("Linearregression",  {
  expect_equal(as.numeric(lr(y ~ x1)[[1]][, 1]), as.numeric(lm(y ~ x1)[[1]]))
  expect_equal(as.numeric(lr(y ~ x1 + x2)[[1]][, 1]), as.numeric(lm(y ~ x1 + x2)[[1]]))
  expect_equal(as.numeric(lr(mtcars$mpg ~ mtcars$cyl + mtcars$disp)[[1]][, 1]), as.numeric(lm(mtcars$mpg ~ mtcars$cyl + mtcars$disp)[[1]]))
  expect_equal(as.numeric(lr(mpg ~ cyl + disp + x4, data = mtcars)[[1]][, 1]), as.numeric(lm(mpg ~ cyl + disp + x4, data = mtcars)[[1]]))
})

##1.2 deliminate intercept in the model
test_that("Linearregression",  {
  expect_equal(as.numeric(lr(y ~ x1 + x2, intercept = FALSE)[[1]][, 1]), as.numeric(lm(y ~ x1 + x2 - 1)[[1]]))
})

##1.3 cope with categorical covariate
test_that("Linearregression",  {
  expect_equal(as.numeric(lr(y ~ x1 + x2 + x3, reference='F')[[1]][, 1]), as.numeric(lm(y ~ x1 + x2 + x3)[[1]]))
})

##1.4 categorical covariate with cell means coding method
test_that("Linearregression",  {
  expect_equal(as.numeric(lr(y ~ x1 + x2 + x3, coding = 'means')[[1]][, 1]), as.numeric(lm(y ~ x1 + x2 + x3 - 1)[[1]]))
})


##1.5 different specific situations under what lr() will give a specified error hint.
##Since test_that cannot test error, we just use comment to show specified situations.
y = c(23, 24, 26, 37, 38, 25, 36, 40)
x5 = c(1, 2, 3, 4, 5, 6, 7, 8)
x6 = c(1, 2, 3, 4, 5, 6, 7, 8)
###1.5.1 lr( ~ x5)
test_that("Linearregression",  {
  expect_error(lr( ~ x5), "The model form is incorrect.", fixed=TRUE)
})
#########with error "The model form is incorrect."
###1.5.2 lr(y ~ x)
test_that("Linearregression",  {
  expect_error(lr(y ~ x), "Cannot find data to fit model.", fixed=TRUE)
})
#########with error "Cannot find data to fit model."
###1.5.3 lr(y ~ x1 + x2 + x3, reference = 'S')
test_that("Linearregression",  {
  expect_error(lr(y ~ x1 + x2 + x3, reference = 'S'), "Reference group was not found.", fixed=TRUE)
})
#########with error "Reference group was not found."
###1.5.4 lr(y ~ x1 + x2 + x3)
test_that("Linearregression",  {
  expect_error(lr(y ~ x1 + x2 + x3), "Multiple predictor 'x's were linear dependent or dimensions of 'x's and y were not matched.", fixed=TRUE)
})
#########with error "Multiple predictor 'x's were linear dependent or dimensions of 'x's and y were not matched."

###2 big samples
set.seed(4)
y1 = rnorm(1e+8,1,10)
x7 = rnorm(1e+8,1,10)
x8 = rnorm(1e+8,1,10)
test_that("Linearregression", {
  expect_equal(as.numeric(lr(y1~x7+x8)[[1]][,1]),as.numeric(lm(y1~x7+x8)[[1]]))
})

