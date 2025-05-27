#
#' @description Default plot preparation method for smooth objects `x'
#' inheriting from "mgcv.smooth".
#' @export
.get_plot_prediction_matrix_and_aux.mgcv.smooth <- function(
    mgcv_term,
    data = NULL,
    n = 100,
    n2 = 40,
    ylim = NULL,
    xlim = NULL,
    too_far = 0.1,
    ...) {
  if (mgcv_term$dim == 1) {
    out <- .get_plot_prediction_matrix_and_aux_plot_smooth_1d(
      mgcv_term = mgcv_term,
      data = data,
      n = n,
      xlim = xlim,
      ...
    )
  }

  if (mgcv_term$dim == 2) {
    out <- .get_plot_prediction_matrix_and_aux_plot_smooth_2d(
      mgcv_term = mgcv_term,
      data = data,
      n2 = n2,
      ylim = ylim,
      xlim = xlim,
      too_far = too_far,
      ...
    )
  }

  if (mgcv_term$dim > 2) {
    out <- .get_plot_prediction_matrix_and_aux_plot_smooth_md(
      mgcv_term = mgcv_term,
      data = data,
      n2 = n2,
      ylim = ylim,
      xlim = xlim,
      too_far = too_far,
      ...
    )
  }

  return(out)
}
