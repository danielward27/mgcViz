#' @description This is the plotting method for random effects (simple random intercepts).
#' @export
get_data.random.effect <- function(term, fitted_terms, trans = identity, ...) {
  P <- prepareP(
    term = term,
    fitted_terms = fitted_terms,
    unconditional = FALSE,
    residuals = TRUE,
    n = 100,
    n2 = NULL,
    lims = NULL,
    too_far = NULL,
    se_with_mean = FALSE
  )
  dat <- list()

  .n <- length(P$fit)
  dat$fit <- data.frame(
    x = qnorm(ppoints(.n)),
    y = sort(P$fit),
    ty = trans(sort(P$fit))
  )

  dat$misc <- list("trans" = trans)
  return(dat)
}


#' @noRd
#' @export
.get_plot_predict_matrix_and_aux.random.effect <- function(
    mgcv_term,
    data = NULL,
    n = 100,
    ...) {
  raw <- data[mgcv_term$term][[1]]
  p <- mgcv_term$last.para - mgcv_term$first.para + 1
  X <- diag(p) # prediction matrix for this term
  list(
    X = X,
    aux = list(raw = raw)
  )
}
