#' @description Get the data for plottingMarkov random field smooths.
#' @export
get_data.mrf.smooth <- function(
    term,
    trans = identity,
    se_with_mean = FALSE,
    unconditional = FALSE,
    ...) {
  p <- prepareP(
    term = term,
    unconditional = unconditional,
    residuals = TRUE,
    res_den = "none",
    n = NULL,
    n2 = NULL,
    ylim = NULL,
    xlim = NULL,
    too_far = NULL,
    se_with_mean = se_with_mean
  )

  data <- list()
  .npb <- sapply(p$smooth$xt$polys, nrow)
  data$fit <- as.data.frame(do.call("rbind", p$smooth$xt$polys))
  names(data$fit) <- c("x", "y")
  data$fit$z <- rep(p$fit, .npb)
  data$fit$tz <- trans(data$fit$z)
  data$fit$id <- rep(p$raw, .npb)
  data$misc <- list("trans" = trans)
  return(data)
}


#' @noRd
#' @export
.get_plot_prediction_matrix_and_aux.mrf.smooth <- function(
    term,
    data,
    ...) {
  raw <- data[term$term][[1]]
  dat <- data.frame(x = factor(names(term$xt$polys), levels = levels(term$knots)))
  names(dat) <- term$term
  X <- PredictMat(term, dat) # prediction matrix for this term
  return(list(
    X = X,
    raw = raw
  ))
}
