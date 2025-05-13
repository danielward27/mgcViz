library(testthat)

test_that("dataframe_extract", {
  library(mgcViz)
  set.seed(2) ## simulate some data...
  dat <- gamSim(1, n = 200, dist = "normal", scale = 2)
  fit <- gam(y ~ s(x0) + s(x1) + s(x2) + s(x3), data = dat)
  p <- plot(sm(getViz(fit), 1)) # TODO currently highjacking plot :D
  expect_true(is.data.frame(p$fit))
  expect_true(is.data.frame(p$res))
})
