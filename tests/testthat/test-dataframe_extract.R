library(testthat)

test_that("1D", {
  library(mgcViz)
  set.seed(2) ## simulate some data...
  dat <- gamSim(1, n = 200)
  fit <- gam(y ~ s(x0), data = dat)
  p <- plot(sm(getViz(fit), 1)) # TODO currently highjacking plot :D
  expect_true(is.data.frame(p$fit))
  expect_true(is.data.frame(p$res))
})


test_that("2D", {
  library(mgcViz)
  set.seed(2) ## simulate some data...
  dat <- gamSim(2, n = 200, verbose = FALSE)$data
  fit <- gam(y ~ s(x, z), data = dat)
  p <- plot(sm(getViz(fit), 1)) # TODO currently highjacking plot :D
  expect_true(is.data.frame(p$fit))
  expect_true(is.data.frame(p$res))
})


test_that("MD", {
  library(mgcViz)
  set.seed(2) ## simulate some data...
  n <- 200
  x <- rnorm(n)
  y <- rnorm(n)
  z <- rnorm(n)
  y <- (x - z)^2 + (y - z)^2 + rnorm(n)
  fit <- gam(y ~ s(x, y, z))
  p <- plot(sm(getViz(fit), 1), fix = c("z" = 0)) # TODO currently highjacking plot :D
  expect_true(is.data.frame(p$fit))
  expect_true(is.data.frame(p$res))
})


test_that("MRF", {
  library(mgcv)
  ## Load Columbus Ohio crime data (see ?columbus for details and credits)
  data(columb) ## data frame
  data(columb.polys) ## district shapes list
  xt <- list(polys = columb.polys) ## neighbourhood structure info for MRF
  fit <- gam(
    crime ~ s(district, bs = "mrf", xt = xt),
    data = columb,
  )
  p <- plot(sm(getViz(fit), 1)) # TODO currently highjacking plot :D
  expect_true(is.data.frame(p$fit))

  # TODO no residuals in mrf. To be expected?
})



test_that("FS Interaction 1D", {
  library(mgcv)
  set.seed(2)
  n <- 200
  group <- factor(rep(letters[1:4], each = 50))
  x <- runif(n)
  y <- sin(2 * pi * x) + as.numeric(group) + rnorm(n, sd = 0.2)
  fit <- gam(
    y ~ s(x, group, bs = "fs", k = 5),
    data = data.frame(y, x, group),
  )
  p <- plot(sm(getViz(fit), 1)) # TODO currently highjacking plot :D
  expect_true(is.data.frame(p$fit))
  # TODO No residuals in FS interaction 1D. To be expected?
})




# TODO just remove trans from all of them?

# What I do
# - Go to private plot
# - Delete plotting code and return data only
# -  Remove class(out) <- c("plotSmooth", "gg")
