#' @rdname get_data.pterm_factor
#' @export get_data.pterm_logical
#' @export
#'
get_data.pterm_logical <- function(term, maxpo = 1e4, trans = identity, ...) {
  get_data.pterm_factor(term = term, maxpo = maxpo, trans = trans, ...)
}
