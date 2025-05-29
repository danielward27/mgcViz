#' @description Get the data for plottingMarkov random field smooths.
#' @export
get_data.mrf.smooth <- function(
    term,
    fitted_terms,
    trans = identity,
    se_with_mean = FALSE,
    unconditional = FALSE,
    ...) {
  p <- prepareP(
    term = term,
    fitted_terms,
    unconditional = unconditional,
    residuals = TRUE,
    n = NULL,
    n2 = NULL,
    lims = NULL,
    too_far = NULL,
    se_with_mean = se_with_mean
  )

  data <- list()
  .npb <- sapply(p$aux$smooth$xt$polys, nrow)
  data$fit <- as.data.frame(do.call("rbind", p$aux$smooth$xt$polys))
  names(data$fit) <- c("x", "y")
  data$fit$z <- rep(p$fit, .npb)
  data$fit$tz <- trans(data$fit$z)
  data$fit$id <- rep(p$aux$raw, .npb)
  data$misc <- list("trans" = trans)
  return(data)
}
#' @noRd
#' @export
.get_plot_predict_matrix_and_aux.mrf.smooth <- function(
    mgcv_term,
    data,
    ...) {
  raw <- data[mgcv_term$term][[1]]
  dat <- data.frame(x = factor(names(mgcv_term$xt$polys), levels = levels(mgcv_term$knots)))
  names(dat) <- mgcv_term$term
  X <- PredictMat(mgcv_term, dat) # prediction matrix for this term
  list(
    X = X,
    aux = list(raw = raw, smooth = mgcv_term)
  )
}
