# Note lims is c(lower, upper) for 1D and list(c(lower1, upper1), c(lower2, upper2)) in 2D/MD cases
#' @export
prepareP <- function(
    term,
    fitted_terms,
    gam,
    unconditional,
    n,
    lims,
    se_with_mean,
    nsim = 0,
    ...) {
  working_residuals <- get_working_residuals(gam)

  term_fit <- fitted_terms[
    , number_parametric(gam) + term$term_idx
  ]
  mgcv_term <- gam$smooth[[term$term_idx]]

  pred_matrix_and_aux <- .get_plot_predict_matrix_and_aux(
    mgcv_term = mgcv_term,
    data = gam$model,
    n = n,
    lims = lims,
    ...
  )

  fit <- .get_fit_plot_data(
    pred_matrix = pred_matrix_and_aux$X,
    term,
    gam
  )

  se_fit <- .get_errors_plot_data(
    gam = gam,
    mgcv_term = mgcv_term,
    pred_matrix = pred_matrix_and_aux$X,
    se_with_mean
  )

  partial_resids <- term_fit + working_residuals

  list(
    fit = fit,
    se_fit = se_fit,
    partial_resids = partial_resids,
    aux = pred_matrix_and_aux$aux
  )
}
