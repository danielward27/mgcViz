# TODO Why doesn't this take the output of perpareP?

#' @export
get_data.fs.interaction.1D <- function(
    term,
    fitted_terms,
    gam = gam,
    n = 100,
    lims = NULL,
    trans = identity,
    ...) {
  P <- prepareP(
    term = term,
    fitted_terms = fitted_terms,
    gam = gam,
    unconditional = FALSE,
    residuals = FALSE,
    n = n,
    n2 = NULL,
    lims = lims,
    se_with_mean = FALSE
  )

  data <- list()
  data$fit <- data.frame(
    "x" = rep(P$aux$x, P$aux$nf),
    "y" = P$fit,
    "ty" = trans(P$fit),
    "id" = as.factor(rep(P$aux$smooth$flev, each = P$aux$n))
  )
  data$misc <- list("trans" = trans)
  data
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
  raw <- data[mgcv_term$base$term][[1]]

  # Generate x sequence for prediction
  if (is.null(lims)) {
    lims <- range(raw)
  }
  x_seq <- seq(lims[1], lims[2], length = n)

  nf <- length(mgcv_term$flev)
  fac <- rep(mgcv_term$flev, rep(n, nf))
  dat <- data.frame(as.factor(fac), x_seq)
  names(dat) <- c(mgcv_term$fterm, mgcv_term$base$term)
  X <- PredictMat(mgcv_term, dat)
  list(
    X = X,
    aux = list(x = x_seq, n = n, nf = nf, smooth = mgcv_term)
  )
}
