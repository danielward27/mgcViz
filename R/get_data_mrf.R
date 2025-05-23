#'
#' Plotting Markov random field smooths
#'
#' @description This is the plotting method for Markov random field smooths.
#' @param x a smooth effect object, extracted using [mgcViz::sm].
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
#' @name plot.mrf.smooth
#' @examples
#' library(mgcViz)
#' ## Load Columbus Ohio crime data (see ?columbus for details and credits)
#' data(columb) ## data frame
#' data(columb.polys) ## district shapes list
#' xt <- list(polys = columb.polys) ## neighbourhood structure info for MRF
#' par(mfrow = c(2, 2))
#' ## First a full rank MRF...
#' b <- gam(crime ~ s(district, bs = "mrf", xt = xt), data = columb, method = "REML")
#' b <- getViz(b)
#'
#' # Manual plot
#' plot(sm(b, 1)) + l_poly(colour = 2) +
#'   scale_fill_gradientn(colours = heat.colors(50))
#'
#' # Default plot
#' plot(b)
#'
#' @rdname plot.mrf.smooth
#' @export get_data.mrf.smooth
#' @export
#'
get_data.mrf.smooth <- function(
    term,
    trans = identity,
    se_with_mean = FALSE,
    unconditional = FALSE,
    ...) {
  p <- prepareP(
    term = term,
    unconditional = unconditional,
    residuals = TRUE,
    resDen = "none",
    se = TRUE,
    se.mult = 1,
    n = NULL,
    n2 = NULL,
    ylim = NULL,
    xlim = NULL,
    too.far = NULL,
    seWithMean = se_with_mean
  )

  data <- list()
  .npb <- sapply(p$smooth$xt$polys, nrow)
  data$fit <- as.data.frame(do.call("rbind", p$smooth$xt$polys))
  names(data$fit) <- c("x", "y")
  data$fit$z <- rep(p$fit, .npb)
  data$fit$tz <- trans(data$fit$z)
  data$fit$id <- rep(p$raw, .npb)
  data$misc <- list("trans" = trans)
  return(data)
}


#' @noRd
#' @export
.prepare.mrf.smooth <- function(
    term,
    data,
    se1.mult,
    se2.mult,
    partial.resids,
    se,
    n,
    n2,
    ylim,
    xlim,
    too.far,
    trans,
    phi,
    theta,
    scheme,
    ...) {
  raw <- data[term$term][[1]]
  dat <- data.frame(x = factor(names(term$xt$polys), levels = levels(term$knots)))
  names(dat) <- term$term
  X <- PredictMat(term, dat) # prediction matrix for this term
  return(list(
    X = X,
    scale = FALSE,
    se = FALSE,
    raw = raw
  ))
}
