# TODO Why doesn't this take the output of perpareP?

#' @export
get_data.fs.interaction.1D <- function(
    term,
    fitted_terms,
    gam,
    n = 100,
    lims = NULL,
    ...) {
  P <- prepareP(
    term = term,
    fitted_terms = fitted_terms,
    gam = gam,
    unconditional = FALSE,
    n = n,
    lims = lims,
    se_with_mean = FALSE
  )
  mgcv_term = gam$smooth[[term$term_idx]]
  fit <- data.frame(
    "x" = rep(P$aux$x, length(mgcv_term$flev)),
    "y" = P$fit,
    "id" = as.factor(rep(mgcv_term$flev, each = length(P$aux$x)))
  )
  list(
    fit = fit
  )
}

#' @noRd
#' @export
.get_plot_predict_matrix_and_aux.fs.interaction <- function(
    mgcv_term,
    data = NULL,
    n = 100,
    lims = NULL,
    ...) {
  if (mgcv_term$dim > 1) {
    stop("no method for base smooth dim > 1")
  }

  # Generate x sequence for prediction
  x_seq <- seq(lims[1], lims[2], length = n)
  num_levels <- length(mgcv_term$flev)
  fac <- rep(mgcv_term$flev, rep(n, num_levels))
  dat <- data.frame(as.factor(fac), x_seq)
  names(dat) <- c(mgcv_term$fterm, mgcv_term$base$term)
  X <- PredictMat(mgcv_term, dat)
  list(
    X = X,
    aux = list(x = x_seq)
  )
}
