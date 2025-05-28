# Prepares P list, smooth (sm) and
#' @export
prepareP <- function(
    term,
    unconditional,
    residuals,
    res_den,
    n,
    n2,
    ylim,
    xlim,
    too_far,
    se_with_mean,
    nsim = 0,
    ...) {
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

  fit_and_errors <- .get_fit_and_errors_dataset(
    term = term,
    unconditional = unconditional,
    residuals = residuals,
    res_den = res_den
  )

  # TODO the appending to the list is not pretty in the function below :)
  plot_data <- .get_fit_and_errors_plot_data(
    pred_matrix_and_aux = pred_matrix_and_aux,
    term = term,
    partial_resids = fit_and_errors$partial_resids,
    se = TRUE, # For now hardcoded but easy to change if needed
    n = n,
    n2 = n2,
    ylim = ylim,
    xlim = xlim,
    too_far = too_far,
    se_with_mean = se_with_mean,
    fit_smooth = fit_and_errors$fv_terms,
    w_resid = fit_and_errors$w_resid,
    res_den = res_den,
    nsim = nsim,
    ...
  )

  list(
    aux = pred_matrix_and_aux$aux,
    fit = plot_data$fit,
    se.fit = plot_data$se.fit,
    p.resid = plot_data$p.resid,
    se = plot_data$se
  )
}

# TODO maybe index and gamviz object is nicer than passing term which implicitly contains gamviz?
