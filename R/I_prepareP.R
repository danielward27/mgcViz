# Prepares P list, smooth (sm) and
#' @export
prepareP <- function(
  term,
  unconditional,
  residuals,
  resDen,
  se,
  se_mult,
  n,
  n2,
  ylim,
  xlim,
  too_far,
  seWithMean,
  nsim = 0,
  ...
) {
  fit_and_errors <- .data_fit_and_errors(
    term = term,
    unconditional = unconditional,
    residuals = residuals,
    resDen = resDen,
    se = se
  )

  plot_data <- .get_plot_data(
    term = term$gam_viz$smooth[[term$term_idx]],
    gam_viz = term$gam_viz,
    partial_resids = fit_and_errors$partial_resids,
    se = fit_and_errors$se,
    n = n,
    n2 = n2,
    ylim = ylim,
    xlim = xlim,
    too_far = too_far,
    se1_mult = se_mult,
    se2_mult = se_mult,
    seWithMean = seWithMean,
    fitSmooth = fit_and_errors$fv_terms,
    w_resid = fit_and_errors$w_resid,
    resDen = resDen,
    nsim = nsim,
    ...
  )

  return(plot_data)  # TODO, I don't like the concat in plot data, we can concat after the fact?
}
