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
    out <- .preparePlotSmooth1D(
      term = term,
      data = data,
      se_mult = se1_mult,
      n = n,
      xlim = xlim,
      ...
    )
  }

  if (term$dim == 2) {
    out <- .preparePlotSmooth2D(
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
    out <- .preparePlotSmoothMD(
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
