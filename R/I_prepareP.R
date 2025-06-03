# Note lims is c(lower, upper) for 1D and list(c(lower1, upper1), c(lower2, upper2)) in 2D/MD cases
#' @export
prepareP <- function(
    term,
    fitted_terms,
    gam,
    unconditional,
    residuals,
    n,
    n2,
    lims,
    se_with_mean,
    nsim = 0,
    ...) {
  working_residuals <- get_working_residuals(gam)

  term_fit <- fitted_terms[
    , number_parametric(gam) + term$term_idx
  ]

  pred_matrix_and_aux <- .get_plot_predict_matrix_and_aux(
    mgcv_term = gam$smooth[[term$term_idx]],
    data = gam$model,
    n = n,
    n2 = n2,
    lims = lims,
    ...
  )

  fit <- .get_fit_plot_data(
    pred_matrix_and_aux,
    term,
    gam
  )

  se_fit <- .get_errors_plot_data(
    gam = gam,
    mgcv_term = gam$smooth[[term$term_idx]],
    pred_and_aux = pred_matrix_and_aux,
    se_with_mean
  )

  partial_resids <- term_fit + working_residuals
  x_raw <- gam$model[gam$smooth[[term$term_idx]]$term][[1]]
  list(
    fit = fit,
    se_fit = se_fit,
    partial_resids = partial_resids,
    aux = pred_matrix_and_aux$aux,
    x_raw = x_raw
  )
}
