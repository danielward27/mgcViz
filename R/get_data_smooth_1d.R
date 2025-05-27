#' @description Get data for plotting one dimensional smooth effects.
#' @export
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
    se_mult = 1,
    n = n,
    n2 = NULL,
    ylim = NULL,
    xlim = xlim,
    too_far = NULL,
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
.prepare_plot_smooth_1d <- function(
    term,
    data,
    se_mult = 1,
    n = 100,
    xlim = NULL,
    ...) {
  out <- NULL
  if (term$plot.me) {
    raw <- as.vector(data[term$term][[1]])
    if (is.null(xlim)) {
      # Generate x sequence for prediction
      x_seq <- seq(min(raw), max(raw), length = n)
    } else {
      x_seq <- seq(xlim[1], xlim[2], length = n)
    }
    if (term$by != "NA") {
      # Deal with any by variables
      by <- rep(1, n)
      dat <- data.frame(x = x_seq, by = by)
      names(dat) <- c(term$term, term$by)
    } else {
      dat <- data.frame(x = x_seq)
      names(dat) <- term$term
    } # Finished preparing prediction data.frame
    X <- PredictMat(term, dat) # prediction matrix for this term

    if (is.null(xlim)) {
      xlim <- range(x_seq)
    }
    out <- list(
      X = X,
      x = x_seq,
      scale = TRUE,
      se = TRUE,
      raw = raw,
      se_mult = se_mult,
      xlim = xlim
    )
  }
  return(out)
}


##
## Default plot preparation method for smooth objects `x' inheriting from "mgcv.smooth"
## Input:
## `x' is a smooth object, usually part of a `gam' fit. It has an attribute
##     'coefficients' containg the coefs for the smooth, but usually these
##     are not needed.
## Output is a list of plot data containing:
##     * fit - the values for plotting
##     * se.fit - standard errors of fit (can be NULL)
##     * the values against which to plot
##     * any raw data information
##     * any partial.residuals

#' @description Default plot preparation method for smooth objects `x'
#' inheriting from "mgcv.smooth".
#' @param x Is a smooth object, usually part of a `gam' fit. It has an attribute
#' `coefficients` containg the coefs for the smooth, but usually these
#' are not needed.
#' @param data The data.
#' @param n Number of points used for each 1-d plot - for a nice smooth plot
#'   this needs to be several times the estimated degrees of freedom for the
#'   smooth. Default value 100.
#' @param n2 Square root of number of points used to grid estimates of 2-d
#'   functions for contouring.
#' @param ylim If supplied then this pair of numbers are used as the y limits
#'   for each plot.
#' @param xlim If supplied then this pair of numbers are used as the x limits
#'   for each plot.
#' @param too_far If greater than 0 then this is used to determine when a
#'   location is too far from data to be plotted when plotting 2-D smooths. This
#'   is useful since smooths tend to go wild away from data. The data are scaled
#'   into the unit square before deciding what to exclude, and too_far <- <-  is a
#'   distance within the unit square. Setting to zero can make plotting faster
#'   for large datasets, but care then needed with interpretation of plots.
#' @param ... Other graphics parameters to pass on to plotting commands.
#' See details for smooth plot specific options.
#' @noRd
#' @export
.prepare.mgcv.smooth <- function(
    term,
    data = NULL,
    se1_mult = 1,
    se2_mult = 2,
    n = 100,
    n2 = 40,
    ylim = NULL,
    xlim = NULL,
    too_far = 0.1,
    ...) {
  if (term$dim == 1) {
    out <- .prepare_plot_smooth_1d(
      term = term,
      data = data,
      se_mult = se1_mult,
      n = n,
      xlim = xlim,
      ...
    )
  }

  if (term$dim == 2) {
    out <- .prepare_plot_smooth_2d(
      term = term,
      data = data,
      se_mult = se2_mult,
      n2 = n2,
      ylim = ylim,
      xlim = xlim,
      too_far = too_far,
      ...
    )
  }

  if (term$dim > 2) {
    out <- .prepare_plot_smooth_md(
      term = term,
      data = data,
      se_mult = se2_mult,
      n2 = n2,
      ylim = ylim,
      xlim = xlim,
      too_far = too_far,
      ...
    )
  }

  return(out)
}
