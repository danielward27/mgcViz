#
#' @description Default plot preparation method for smooth objects `x'
#' inheriting from "mgcv.smooth".
#' @export
.get_plot_prediction_matrix_and_aux.mgcv.smooth <- function(
    term,
    data = NULL,
    se1_mult = 1,
    se2_mult = 2,
    n = 100,
    n2 = 40,
    ylim = NULL,
    xlim = NULL,
    too_far = 0.1,
    ...) {
  if (term$dim == 1) {
    out <- .get_plot_prediction_matrix_and_aux_plot_smooth_1d(
      term = term,
      data = data,
      se_mult = se1_mult,
      n = n,
      xlim = xlim,
      ...
    )
  }

  if (term$dim == 2) {
    out <- .get_plot_prediction_matrix_and_aux_plot_smooth_2d(
      term = term,
      data = data,
      se_mult = se2_mult,
      n2 = n2,
      ylim = ylim,
      xlim = xlim,
      too_far = too_far,
      ...
    )
  }

  if (term$dim > 2) {
    out <- .get_plot_prediction_matrix_and_aux_plot_smooth_md(
      term = term,
      data = data,
      se_mult = se2_mult,
      n2 = n2,
      ylim = ylim,
      xlim = xlim,
      too_far = too_far,
      ...
    )
  }

  return(out)
}
