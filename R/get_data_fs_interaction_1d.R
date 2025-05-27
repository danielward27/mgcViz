# TODO Why doesn't this take the output of perpareP?

#' @export
get_data.fs.interaction.1D <- function(
    term,
    n = 100,
    xlim = NULL,
    trans = identity,
    ...) {
  P <- prepareP(
    term = term,
    unconditional = FALSE,
    residuals = FALSE,
    res_den = "none",
    n = n,
    n2 = NULL,
    ylim = NULL,
    xlim = xlim,
    too_far = NULL,
    se_with_mean = FALSE
  )

  data <- list()
  data$fit <- data.frame(
    "x" = rep(P$x, P$nf),
    "y" = P$fit,
    "ty" = trans(P$fit),
    "id" = as.factor(rep(P$smooth$flev, each = P$n))
  )
  data$misc <- list("trans" = trans)
  return(data)
}



#' @noRd
#' @export
.get_plot_prediction_matrix_and_aux.fs.interaction <- function(
    term,
    data = NULL,
    n = 100,
    ylim = NULL,
    xlim = NULL,
    ...) {
  if (term$dim > 1) {
    stop("no method for base smooth dim > 1")
  }
  raw <- data[term$base$term][[1]]

  # Generate x sequence for prediction
  if (is.null(xlim)) {
    xlim <- range(raw)
  }
  x_seq <- seq(xlim[1], xlim[2], length = n)

  nf <- length(term$flev)
  fac <- rep(term$flev, rep(n, nf))
  dat <- data.frame(as.factor(fac), x_seq)
  names(dat) <- c(term$fterm, term$base$term)
  X <- PredictMat(term, dat)
  return(list(
    X = X,
    raw = raw,
    x = x_seq,
    n = n,
    nf = nf
  ))
}
