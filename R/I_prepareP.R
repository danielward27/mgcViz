#' @export
prepareP <- function(
    term,
    unconditional,
    residuals,
    n,
    n2,
    ylim,
    xlim,
    too_far,
    se_with_mean,
    nsim = 0,
    ...) {
  covariance_and_errors <- .get_covariance_and_errors(
    term = term,
    unconditional = unconditional
  )

  # TODO unused...
  term_fit <- term$gam_viz$store$termsFit[
    , term$gam_viz$store$np + term$term_idx
  ]

  pred_matrix_and_aux <- .get_plot_predict_matrix_and_aux(
    mgcv_term = term$gam_viz$smooth[[term$term_idx]],
    data = term$gam_viz$model,
    n = n,
    n2 = n2,
    ylim = ylim,
    xlim = xlim,
    too_far = too_far,
    ...
  )

  plot_data <- .get_fit_and_errors_plot_data(
    pred_matrix_and_aux = pred_matrix_and_aux,
    term = term,
    compute_partial_resids = TRUE, # For now hardcoded
    compute_se = TRUE, # For now hardcoded
    se_with_mean = se_with_mean,
    term_fit = term_fit, # TODO inconsistent naming
    w_resid = covariance_and_errors$w_resid,
    nsim = nsim
  )

  list(
    fit = plot_data$fit,
    se.fit = plot_data$se.fit,
    partial_resids = plot_data$partial_resids,
    se = plot_data$se,
    aux = pred_matrix_and_aux$aux
  )
}

# TODO maybe index and gamviz object is nicer than passing term which implicitly contains gamviz?
