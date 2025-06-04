#' @description Get the data for plottingMarkov random field smooths.
#' @export
get_data.mrf.smooth <- function(
    term,
    fitted_terms,
    gam = gam,
    se_with_mean = FALSE,
    unconditional = FALSE,
    ...) {
  p <- prepareP(
    term = term,
    fitted_terms,
    gam = gam,
    unconditional = unconditional,
    n = NULL,
    lims = NULL,
    se_with_mean = se_with_mean
  )
  smooth <- gam$smooth[[term$term_idx]]
  .npb <- sapply(smooth$xt$polys, nrow)
  fit <- as.data.frame(do.call("rbind", smooth$xt$polys))
  names(fit) <- c("x1", "x2")
  fit$y <- rep(p$fit, .npb)
  fit$id <- rep(p$x_raw, .npb)
  list(
    fit = fit
  )
}

#' @noRd
#' @export
.get_plot_predict_matrix_and_aux.mrf.smooth <- function(
    mgcv_term,
    data,
    ...) {
  dat <- data.frame(x = factor(names(mgcv_term$xt$polys), levels = levels(mgcv_term$knots)))
  names(dat) <- mgcv_term$term
  X <- PredictMat(mgcv_term, dat) # prediction matrix for this term
  list(
    X = X,
    aux = list()
  )
}
