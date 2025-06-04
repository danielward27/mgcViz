#' @description Get data for plotting smooth effects on the sphere.

#' @export
#'
get_data.sos.smooth <- function(
    term,
    fitted_terms,
    gam = gam,
    n = 40,
    lims = NULL,
    se_with_mean = FALSE,
    unconditional = FALSE,
    ...) {

  P <- prepareP(
    term = term,
    fitted_terms = fitted_terms,
    gam = gam,
    unconditional = unconditional,
    n = n,
    lims = lims,
    se_with_mean = se_with_mean,
  )

  .get_data_shared_2d(
    P = P,
    flip = TRUE
  )
}

#' @noRd
#' @export
.get_plot_predict_matrix_and_aux.sos.smooth <- function(
    mgcv_term,
    data,
    n,
    lims,
    ...) {
    .get_plot_predict_matrix_and_aux_plot_smooth_2d(
      mgcv_term = mgcv_term,
      data = data,
      n = n,
      lims = lims
    )
  }
