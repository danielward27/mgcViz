#'
#' Plotting two dimensional smooth effects
#'
#' @description Plotting method for two dimensional smooth effects.
#'
#' @param x a smooth effect object, extracted using [mgcViz::sm].
#' @param n sqrt of the number of grid points used to compute the effect plot.
#' @param xlim if supplied then this pair of numbers are used as the x limits for the plot.
#' @param ylim if supplied then this pair of numbers are used as the y limits for the plot.
#' @param maxpo maximum number of residuals points that will be used by layers such as
#'              \code{resRug()} and \code{resPoints()}. If number of datapoints > \code{maxpo},
#'              then a subsample of \code{maxpo} points will be taken.
#' @param too.far if greater than 0 then this is used to determine when a location is too far
#'               from data to be plotted. This is useful since smooths tend to go wild
#'               away from data. The data are scaled into the unit square before deciding
#'               what to exclude, and too.far is a distance within the unit square.
#'               Setting to zero can make plotting faster for large datasets, but care
#'               then needed with interpretation of plots.
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
#' @param ... currently unused.
#' @return An objects of class \code{plotSmooth}.
#' @references Marra, G and S.N. Wood (2012) Coverage Properties of Confidence Intervals for
#'             Generalized Additive Model Components. Scandinavian Journal of Statistics.
#' @name plot.mgcv.smooth.2D
#' @examples
#' library(mgcViz)
#' set.seed(2) ## simulate some data...
#' dat <- gamSim(1, n = 700, dist = "normal", scale = 2)
#' b <- gam(y ~ s(x0) + s(x1, x2) + s(x3), data = dat, method = "REML")
#' b <- getViz(b)
#'
#' # Plot 2D effect with noised-up raster, contour and rug for design points
#' # Opacity is proportional to the significance of the effect
#' plot(sm(b, 2)) + l_fitRaster(pTrans = zto1(0.05, 2, 0.1), noiseup = TRUE) +
#'   l_rug() + l_fitContour()
#'
#' # Plot contour of effect joint density of design points
#' plot(sm(b, 2)) + l_dens(type = "joint") + l_points() + l_fitContour() +
#'   coord_cartesian(expand = FALSE) # Fill the plot
#'
#' ###
#' # Quantile GAM example
#' ###
#' b <- mqgamV(y ~ s(x0) + s(x1, x2) + s(x3), qu = c(0.3, 0.7), data = dat)
#'
#' plot(sm(b, 2)) + l_fitRaster(noiseup = TRUE) + l_fitContour(colour = 2)
#'
#' @importFrom mgcv exclude.too.far
#' @rdname plot.mgcv.smooth.2D
#' @export get_data.mgcv.smooth.2D
#' @export
#'
get_data.mgcv.smooth.2D <- function(
    term,
    n = 40,
    xlim = NULL,
    ylim = NULL,
    maxpo = 1e4,
    too.far = 0.1,
    trans = identity,
    seWithMean = FALSE,
    unconditional = FALSE,
    ...) {
  # 1) Prepare data
  P <- prepareP(
    term = term,
    unconditional = unconditional,
    residuals = TRUE,
    resDen = "none",
    se = TRUE,
    se.mult = 1,
    n = NULL,
    n2 = n,
    ylim = ylim,
    xlim = xlim,
    too.far = too.far,
    seWithMean = seWithMean
  )

  # 2) Produce output object
  out <- .get_data_shared_2d(term = P$smooth, P = P, trans = trans, maxpo = maxpo)
  return(out)
}


# Used by e.g MD and SOS
#' @noRd
.get_data_shared_2d <- function(term, P, trans, maxpo, flip = FALSE) {  # TODO is term even used?
  .dat <- list()
  # 1) Build dataset on fitted effect
  P$fit[P$exclude] <- NA
  .dat$fit <- data.frame(
    "z" = drop(P$fit),
    "tz" = drop(trans(P$fit)),
    "x" = rep(P$x, length(P$fit) / length(P$x)),
    "y" = rep(P$y, each = length(P$fit) / length(P$x)),
    "se" = P$se
  )

  # 2) Build dataset on residuals
  if (!is.null(P$raw)) {
    # Exclude points too far from current slice (relevant only when called by plot.mgcv.smooth.MD)
    if (!is.null(P$exclude2) && any(P$exclude2)) {
      P$raw <- P$raw[!P$exclude2, ]
    }

    # Exclude residuals falling outside boundaries
    .dat$res <- P$raw[
      P$raw$x >= P$xlim[1] &
        P$raw$x <= P$xlim[2] &
        P$raw$y >= P$ylim[1] &
        P$raw$y <= P$ylim[2], ,
      drop = FALSE
    ]

    # Sample if too many points (> maxpo)
    nres <- nrow(.dat$res)
    .dat$res$sub <- if (nres > maxpo) {
      sample(c(rep(T, maxpo), rep(F, nres - maxpo)))
    } else {
      rep(T, nres)
    }
  }

  .dat$misc <- list("trans" = trans)

  if (flip) {
    .dat2 <- .dat
    .dat2$fit$x <- .dat$fit$y
    .dat2$fit$y <- .dat$fit$x
    .dat2$res$x <- .dat$res$y
    .dat2$res$y <- .dat$res$x
    .dat <- .dat2
  }

  return(.dat)
}


# Internal function for preparing plot of two dimensional smooths
.preparePlotSmooth2D <- function(
    term,
    data = NULL,
    se.mult = 2,
    n2 = 40,
    ylim = NULL,
    xlim = NULL,
    too.far = 0.1,
    ...) {
  out <- NULL
  if (term$plot.me) {
    xterm <- term$term[1]
    yterm <- term$term[2]
    raw <- data.frame(
      x = as.numeric(data[xterm][[1]]),
      y = as.numeric(data[yterm][[1]])
    )
    n2 <- max(10, n2)
    if (is.null(xlim)) {
      xm <- seq(min(raw$x), max(raw$x), length = n2)
    } else {
      xm <- seq(xlim[1], xlim[2], length = n2)
    }
    if (is.null(ylim)) {
      ym <- seq(min(raw$y), max(raw$y), length = n2)
    } else {
      ym <- seq(ylim[1], ylim[2], length = n2)
    }
    xx <- rep(xm, n2)
    yy <- rep(ym, rep(n2, n2))
    if (too.far > 0) {
      exclude <- exclude.too.far(xx, yy, raw$x, raw$y, dist = too.far)
    } else {
      exclude <- rep(FALSE, n2 * n2)
    }
    if (term$by != "NA") {
      # deal with any by variables
      by <- rep(1, n2^2)
      dat <- data.frame(x = xx, y = yy, by = by)
      colnames(dat) <- c(xterm, yterm, term$by)
    } else {
      dat <- data.frame(x = xx, y = yy)
      colnames(dat) <- c(xterm, yterm)
    } ## prediction data.frame complete
    X <- PredictMat(term, dat) ## prediction matrix for this term
    
    if (is.null(ylim)) {
      ylim <- range(ym)
    }
    if (is.null(xlim)) {
      xlim <- range(xm)
    }
    out <- list(
      X = X,
      x = xm,
      y = ym,
      scale = FALSE,
      se = TRUE,
      raw = raw,
      se.mult = se.mult,
      ylim = ylim,
      xlim = xlim,
      exclude = exclude
    )
  }
  return(out)
}
