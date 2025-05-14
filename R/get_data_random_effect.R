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
get_data.random.effect <- function(x, trans = identity, ...) {
  P <- prepareP(
    o = x,
    unconditional = FALSE,
    residuals = TRUE,
    resDen = "none",
    se = TRUE,
    se.mult = 1,
    n = 100,
    n2 = NULL,
    xlab = NULL,
    ylab = NULL,
    main = NULL,
    ylim = NULL,
    xlim = NULL,
    too.far = NULL,
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
    x,
    data = NULL,
    label = "",
    n = 100,
    xlab = NULL,
    ylab = NULL,
    main = NULL,
    ylim = NULL,
    xlim = NULL,
    ...) {
  raw <- data[x$term][[1]]
  p <- x$last.para - x$first.para + 1
  X <- diag(p) # prediction matrix for this term
  if (is.null(xlab)) xlabel <- "Gaussian quantiles" else xlabel <- xlab
  if (is.null(ylab)) ylabel <- "effects" else ylabel <- ylab
  if (!is.null(main)) label <- main

  return(list(
    X = X,
    scale = FALSE,
    se = FALSE,
    raw = raw,
    xlab = xlabel,
    ylab = ylabel,
    main = label
  ))
}
