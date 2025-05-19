#' @rdname get_data.ptermFactor
#' @export get_data.ptermLogical
#' @export
#'
get_data.ptermLogical <- function(term, maxpo = 1e4, trans = identity, ...) {
  get_data.ptermFactor(term = term, maxpo = maxpo, trans = trans, ...)
}
