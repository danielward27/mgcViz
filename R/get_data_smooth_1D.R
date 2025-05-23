#'
#' Plotting one dimensional smooth effects
#'
#' @description Plotting method for one dimensional smooth effects.
#' @name plot.mgcv.smooth.1D
#' @param x a smooth effect object, extracted using [mgcViz::sm].
#' @param n number of grid points used to compute main effect and c.i. lines.
#'          For a nice smooth plot this needs to be several times the estimated degrees of
#'          freedom for the smooth.
#' @param xlim if supplied then this pair of numbers are used as the x limits for the plot.
#' @param maxpo maximum number of residuals points that will be used by layers such as
#'              \code{resRug()} and \code{resPoints()}. If number of datapoints > \code{maxpo},
#'              then a subsample of \code{maxpo} points will be taken.
#' @param trans monotonic function to apply to the smooth and residuals, before plotting.
#'              Monotonicity is not checked.
#' @param unconditional if \code{TRUE} then the smoothing parameter uncertainty corrected covariance
#'                      matrix is used to compute uncertainty bands, if available.
#'                      Otherwise the bands treat the smoothing parameters as fixed.
#' @param seWithMean if TRUE the component smooths are shown with confidence intervals that
#'                   include the uncertainty about the overall mean. If FALSE then the uncertainty
#'                   relates purely to the centred smooth itself. Marra and Wood (2012) suggests
#'                   that TRUE results in better coverage performance, and this is also suggested
#'                   by simulation.
#' @param nsim number of smooth effect curves to be simulated from the posterior distribution.
#'             These can be plotted using the \link{l_simLine} layer. See Examples section below.
#' @param asFact determines whether to use a factor or colour bar legend for plot.multi.mgcv.smooth.1D.
#'               For most models the default is \code{TRUE}. When working with QGAM models fitted with
#'               \link{mqgamV}, the default is \code{FALSE} for less than 10 quantiles,
#'               \code{TRUE} otherwise. For QGAM models there a third option, \code{asFact = "force"},
#'               which forces the use of a discrete colour scale.
#' @param ... currently unused.
#' @return An objects of class \code{plotSmooth}.
#' @references Marra, G and S.N. Wood (2012) Coverage Properties of Confidence Intervals for
#'             Generalized Additive Model Components. Scandinavian Journal of Statistics.
#' @examples
#' library(mgcViz)
#' n <- 1e3
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' dat <- data.frame(
#'   "x1" = x1, "x2" = x2,
#'   "y" = sin(x1) + 0.5 * x2^2 + pmax(x2, 0.2) * rnorm(n)
#' )
#' b <- bamV(y ~ s(x1) + s(x2), data = dat, method = "fREML", aGam = list(discrete = TRUE))
#'
#' o <- plot(sm(b, 1), nsim = 50) # 50 posterior simulations
#'
#' # Plot with fitted effect + posterior simulations + rug on x axis
#' (o <- o + l_simLine() + l_fitLine(colour = "red") +
#'   l_rug(alpha = 0.8))
#'
#' # Add CI lines at 1*sigma and 5*sigma
#' (o <- o + l_ciLine(mul = 1) + l_ciLine(mul = 5, colour = "blue", linetype = 2))
#'
#' # Add partial residuals and change theme
#' (o + l_points(shape = 19, size = 1, alpha = 0.2) + theme_classic())
#'
#' # Get second effect plot
#' o2 <- plot(sm(b, 2))
#'
#' # Plot it with polygon for partial residuals
#' o2 + l_ciPoly(mul = 5, fill = "light blue") +
#'   l_fitLine(linetype = 2, colour = "red")
#'
#' # Plot is with conditional density of partial residuals
#' o2 + l_dens(type = "cond", alpha = 0.9) +
#'   l_fitLine(linetype = 2, colour = "red")
#'
#' ########
#' # Quantile GAM example
#' ########
#' # Fit model
#' b <- mqgamV(y ~ s(x1) + s(x2), qu = c(0.2, 0.5, 0.8), data = dat)
#'
#' plot(sm(b, 1)) + l_fitLine(linetype = 2) + l_rug(colour = "blue")
#'
#' @rdname plot.mgcv.smooth.1D
#' @export get_data.mgcv.smooth.1D
#' @export
#'
get_data.mgcv.smooth.1D <- function(
    term,
    n = 100,
    xlim = NULL,
    maxpo = 1e4,
    trans = identity,
    unconditional = FALSE,
    seWithMean = FALSE,
    nsim = 0,
    ...) {
  P <- prepareP(
    term = term,
    unconditional = unconditional,
    residuals = TRUE,
    resDen = "none",
    se = TRUE,
    se.mult = 1,
    n = n,
    n2 = NULL,
    ylim = NULL,
    xlim = xlim,
    too.far = NULL,
    seWithMean = seWithMean,
    nsim = nsim
  )
  .dat <- list()
  if (!is.null(P$raw)) {
    # Construct data.frame of partial residuals
    res <- data.frame("x" = as.vector(P$raw))
    if (!is.null(P$p.resid) & length(P$p.resid)) {
      res$y <- trans(P$p.resid)
    }

    # Exclude residuals falling outside boundaries
    .dat$res <- res[res$x >= P$xlim[1] & res$x <= P$xlim[2], , drop = FALSE]

    # Sample if too many points (> maxpo)
    nres <- nrow(.dat$res)
    .dat$res$sub <- if (nres > maxpo) {
      sample(c(rep(T, maxpo), rep(F, nres - maxpo)))
    } else {
      rep(T, nres)
    }
  }

  .dat$fit <- data.frame(
    x = P$x, # x values
    y = P$fit, # fitted values
    ty = trans(P$fit), # fitted values after trans
    se = P$se
  ) # standard error

  if (!is.null(P$simF)) {
    nsim <- ncol(P$simF)
    .dat$sim <- data.frame(
      "x" = rep(P$x, nsim),
      "ty" = trans(as.vector(P$simF)),
      "id" = as.factor(rep(1:nsim, each = nrow(P$simF)))
    )
  }

  .dat$misc <- list("trans" = trans)
  return(.dat)
}


# Internal function for preparing plot of one dimensional smooths
.preparePlotSmooth1D <- function(
    term,
    data,
    se.mult = 1,
    n = 100,
    xlim = NULL,
    ...) {
  out <- NULL
  if (term$plot.me) {
    raw <- as.vector(data[term$term][[1]])
    if (is.null(xlim)) {
      # Generate x sequence for prediction
      xx <- seq(min(raw), max(raw), length = n)
    } else {
      xx <- seq(xlim[1], xlim[2], length = n)
    }
    if (term$by != "NA") {
      # Deal with any by variables
      by <- rep(1, n)
      dat <- data.frame(x = xx, by = by)
      names(dat) <- c(term$term, term$by)
    } else {
      dat <- data.frame(x = xx)
      names(dat) <- term$term
    } # Finished preparing prediction data.frame
    X <- PredictMat(term, dat) # prediction matrix for this term

    if (is.null(xlim)) {
      xlim <- range(xx)
    }
    out <- list(
      X = X,
      x = xx,
      scale = TRUE,
      se = TRUE,
      raw = raw,
      se.mult = se.mult,
      xlim = xlim
    )
  }
  return(out)
}
