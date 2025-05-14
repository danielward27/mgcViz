library(testthat)

test_that("1D", {
  library(mgcViz)
  set.seed(2) ## simulate some data...
  dat <- gamSim(1, n = 200, verbose = FALSE)
  fit <- gam(y ~ s(x0), data = dat)
  data <- get_data(sm(getViz(fit), 1))
  expect_true(is.data.frame(data$fit))
  expect_true(is.data.frame(data$res))
})


test_that("2D", {
  library(mgcViz)
  set.seed(2) ## simulate some data...
  dat <- gamSim(2, n = 200, verbose = FALSE)$data
  fit <- gam(y ~ s(x, z), data = dat)
  data <- get_data(sm(getViz(fit), 1))
  expect_true(is.data.frame(data$fit))
  expect_true(is.data.frame(data$res))
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
  data <- get_data(sm(getViz(fit), 1), fix = c("z" = 0))
  expect_true(is.data.frame(data$fit))
  expect_true(is.data.frame(data$res))
})


test_that("MRF", {
  library(mgcv)
  data(columb)
  data(columb.polys)
  xt <- list(polys = columb.polys) ## neighbourhood structure info for MRF
  fit <- gam(
    crime ~ s(district, bs = "mrf", xt = xt),
    data = columb,
  )
  data <- get_data(sm(getViz(fit), 1))
  expect_true(is.data.frame(data$fit))
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
  data <- get_data(sm(getViz(fit), 1))
  expect_true(is.data.frame(data$fit))
  # TODO No residuals in FS interaction 1D. To be expected?
})



test_that("SOS", {
  library(mgcv)
  set.seed(2)
  n <- 200

  f <- function(la, lo) {
    sin(lo) * cos(la - .3)
  }
  lo <- runif(n) * 2 * pi - pi
  la <- runif(3 * n) * pi - pi / 2
  ind <- runif(3 * n) <= cos(la)
  la <- la[ind]
  la <- la[1:n]

  ff <- f(la, lo)
  y <- ff + 0.2 * rnorm(n)

  dat <- data.frame(la = la * 180 / pi, lo = lo * 180 / pi, y = y)
  fit <- gam(y ~ s(la, lo, bs = "sos", k = 60), data = dat)
  data <- get_data(sm(getViz(fit), 1))

  expect_true(is.data.frame(data$fit))
  expect_true(is.data.frame(data$res))
})


test_that("P Term Factor", {
  library(mgcv)
  set.seed(2)
  dat <- gamSim(1, n = 200, dist = "normal", scale = 10, verbose = FALSE)
  dat$fac <- as.factor(sample(c("A1", "A2", "A3"), nrow(dat), replace = TRUE))
  fit <- gam(y ~ fac, data = dat)
  data <- get_data(pterm(getViz(fit), 1)) # Note parametric terms have to use pterm.
  expect_true(is.data.frame(data$fit))
  expect_true(is.data.frame(data$res))
})



test_that("P Term Logical", {
  library(mgcv)
  set.seed(2)
  dat <- gamSim(1, n = 200, dist = "normal", scale = 10, verbose = FALSE)
  dat$logical <- as.logical(sample(c(TRUE, FALSE), nrow(dat), replace = TRUE))
  fit <- gam(y ~ logical, data = dat)
  data <- get_data(pterm(getViz(fit), 1)) # Note parametric terms have to use pterm.
  expect_true(is.data.frame(data$fit))
  expect_true(is.data.frame(data$res))
})



test_that("P Term Numeric", {
  library(mgcv)
  set.seed(2)
  dat <- gamSim(1, n = 200, dist = "normal", scale = 10, verbose = FALSE)
  dat$numeric <- rnorm(200)
  fit <- gam(y ~ numeric, data = dat)
  data <- get_data(pterm(getViz(fit), 1)) # Note parametric terms have to use pterm.
  expect_true(is.data.frame(data$fit))
  expect_true(is.data.frame(data$res))
})


test_that("Random effect", {
  library(mgcViz)
  set.seed(2)
  fit <- gam(travel ~ s(Rail, bs = "re"), data = Rail, method = "REML")

  data <- get_data(sm(getViz(fit), 1))
  expect_true(is.data.frame(data$fit))
  # TODO No residuals for random effects. Expected?
})
