#
#' @description Default plot preparation method for smooth objects `x'
#' inheriting from "mgcv.smooth".
#' @export
.get_plot_predict_matrix_and_aux.mgcv.smooth <- function(
    mgcv_term,
    data = NULL,
    n = 100,
    lims = NULL,
    ...) {
  if (mgcv_term$dim == 1) {
    out <- .get_plot_predict_matrix_and_aux_plot_smooth_1d(
      mgcv_term = mgcv_term,
      data = data,
      n = n,
      lims = lims,
      ...
    )
  }

  if (mgcv_term$dim == 2) {
    out <- .get_plot_predict_matrix_and_aux_plot_smooth_2d(
      mgcv_term = mgcv_term,
      data = data,
      n = n,
      lims = lims,
      ...
    )
  }

  if (mgcv_term$dim > 2) {
    out <- .get_plot_predict_matrix_and_aux_plot_smooth_md(
      mgcv_term = mgcv_term,
      data = data,
      n = n,
      lims = lims,
      ...
    )
  }
  out
}
