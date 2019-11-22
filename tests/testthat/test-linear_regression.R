test_that("multiplication works", {
  expect_equal(2 * 2, 4)
  y=c(1,2,3,4,5)
  x1=c(1,2,3,4,5)
  x2=c(3,4,2,5,2)
  expect_equal(linear_regression(y~x1+x2),2)
})
