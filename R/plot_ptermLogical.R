#' @rdname plot.ptermFactor
#' @export get_data.ptermLogical
#' @export
#'
plot.ptermLogical <- function(x, maxpo = 1e4, trans = identity, ...) {
  plot.ptermFactor(x = x, maxpo = maxpo, trans = trans, ...)
}
