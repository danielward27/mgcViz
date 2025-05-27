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
  fit_and_errors <- .get_fit_and_errors_dataset(
    term = term,
    unconditional = unconditional,
    residuals = residuals,
    res_den = res_den
  )



  plot_data <- .get_fit_and_errors_plot_data(
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

  plot_data # TODO, I don't like the concat in plot data, we can concat after the fact?
}
