#' @description Get the data for plotting of boolean effects
#' @export
#'
get_data.pterm_logical <- function(term, trans = identity, ...) {
  get_data.pterm_factor(term = term, trans = trans, ...)
}
