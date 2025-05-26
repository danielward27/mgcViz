#'
#' Plotting random effects
#'
#' @description This is the plotting method for random effects (simple random intercepts).
#' @param x a random effect object, extracted using [mgcViz::sm].
#' @param trans monotonic function to apply to the fit, confidence intervals and residuals,
#'              before plotting. Monotonicity is not checked.
#' @param ... currently unused.
#' @return An object of class \code{plotSmooth}.
#' @name plot.random.effect
#' @rdname plot.random.effect
#' @export get_data.random.effect
#' @export
#' @examples
#' library(mgcViz)
#' b <- gam(travel ~ s(Rail, bs = "re"), data = Rail, method = "REML")
#' b <- getViz(b)
#' plot(sm(b, 1)) + l_fitLine(colour = 2, linetype = 2) + l_points() +
#'   l_ciLine(colour = 4, linetype = 3)
#'
#' plot(sm(b, 1)) + l_ciPoly() + l_points()
#'
#' # Default
#' plot(b)
#'
#' ###
#' # Quantile GAM version
#' ###
#' b <- mqgamV(travel ~ s(Rail, bs = "re"), data = as.data.frame(Rail), qu = c(0.2, 0.4, 0.6, 0.8))
#'
#' plot(sm(b, 1)) + l_ciPoly() + l_points()
#'
#' # Default
#' plot(b)
#'
get_data.random.effect <- function(term, trans = identity, ...) {
  P <- prepareP(
    term = term,
    unconditional = FALSE,
    residuals = TRUE,
    resDen = "none",
    se = TRUE,
    se_mult = 1,
    n = 100,
    n2 = NULL,
    ylim = NULL,
    xlim = NULL,
    too_far = NULL,
    seWithMean = FALSE
  )
  dat <- list()

  .n <- length(P$fit)
  dat$fit <- data.frame(
    x = qnorm(ppoints(.n)),
    y = sort(P$fit),
    ty = trans(sort(P$fit))
  )

  dat$misc <- list("trans" = trans)
  return(dat)
}


#' @noRd
#' @export
.prepare.random.effect <- function(
    term,
    data = NULL,
    n = 100,
    ylim = NULL,
    xlim = NULL,
    ...) {
  raw <- data[term$term][[1]]
  p <- term$last.para - term$first.para + 1
  X <- diag(p) # prediction matrix for this term

  return(list(
    X = X,
    scale = FALSE,
    se = FALSE,
    raw = raw
  ))
}
