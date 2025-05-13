#'
#' Generic plotting of differences
#'
#' @description Generic function for plotting differences between objects. Useful
#'              mainly for plotting the differences between two smooth effects.
#' @param ... arguments to be passed to methods. This first one will determine which
#'            method will be called.
#' @seealso [plotDiff.mgcv.smooth.1D], [plotDiff.mgcv.smooth.2D], [plotDiff.sos.smooth]
#' @rdname plotDiff
#' @export plotDiff
plotDiff <- function(...) UseMethod("plotDiff")

#'
#' Generic QQ plots
#'
#' @description Generic function for producing QQ-plots.
#' @param ... arguments to be passed to methods. This first one will determine which
#'            method will be called.
#' @seealso [qq.gamViz]
#' @rdname qq
#' @export qq
qq <- function(...) UseMethod("qq")

###### Internal generics
.prepare <- function(...) UseMethod(".prepare")
