# library(testthat)

# test_that("Test order of pterm is consistent with model.", {
#     library(mgcViz)
#     set.seed(2) ## simulate some data...


#     dat <- data.frame(
#         factor = factor(sample(c("A", "B", "C"), size = 200, replace = TRUE)),
#         bool = sample(c(TRUE, FALSE), size = 200, replace = TRUE),
#         x = runif(200, min = 0, max = 100),
#         y = runif(200, min = 0, max = 100)
#     )


#     fit <- gam(y ~ x + s(x) + x:y + bool + factor, data = dat)

#     viz <- getViz(fit)

#     data <- get_data(sm(, 1))


#     expect_true(is.data.frame(data$fit))
#     expect_true(is.data.frame(data$res))
# })
