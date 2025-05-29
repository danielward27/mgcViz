#' @description Get the data for plotting of boolean effects
#' @export
#'
get_data.pterm_logical <- function(term, fitted_terms, trans = identity, ...) {
  get_data.pterm_factor(
    term = term,
    fitted_terms = fitted_terms,
    trans = trans,
    ...
  )
}
