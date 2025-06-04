#' @description This is the plotting method for random effects (simple random intercepts).
#' @export
get_data.random.effect <- function(term, fitted_terms, gam, ...) {
  P <- prepareP(
    term = term,
    fitted_terms = fitted_terms,
    gam = gam,
    unconditional = FALSE,
    n = 100,
    lims = NULL,
    se_with_mean = FALSE
  )
  dat <- list()

  .n <- length(P$fit)
  dat$fit <- data.frame(
    x = qnorm(ppoints(.n)),
    y = sort(P$fit)
  )
  dat
}


#' @noRd
#' @export
.get_plot_predict_matrix_and_aux.random.effect <- function(
    mgcv_term,
    data = NULL,
    n = 100,
    ...) {
  p <- mgcv_term$last.para - mgcv_term$first.para + 1
  X <- diag(p) # prediction matrix for this term
  list(
    X = X,
    aux = list()
  )
}
