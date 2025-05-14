#' @rdname get_data.ptermFactor
#' @export get_data.ptermLogical
#' @export
#'
get_data.ptermLogical <- function(x, maxpo = 1e4, trans = identity, ...) {
  get_data.ptermFactor(x = x, maxpo = maxpo, trans = trans, ...)
}
