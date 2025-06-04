#' @description Get the data for plotting of boolean effects
#' @export
#'
get_data.pterm_logical <- function(term, fitted_terms, ...) {
  get_data.pterm_factor(
    term = term,
    fitted_terms = fitted_terms,
    ...
  )
}
