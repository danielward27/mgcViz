.get_covariance_and_errors <- function(
    term,
    unconditional) { # TODO still doing two relatively unrelated actions...
  # Use Bayesian cov matrix including smoothing parameter uncertainty?
  covariance <- NULL
  if (unconditional) {
    if (is.null(term$gam$Vc)) {
      warning("Smoothness uncertainty corrected covariance not available")
    } else {
      covariance <- term$gam$Vc
    }
  }

  # TODO more descriptive name for residuals, looks like below we get Pearson residuals,
  # not absolute.

  # produce working residuals if info available
  if (is.null(term$gam$residuals) || is.null(term$gam$weights)) {
    stop("Cannot compute working residuals.")
  }
  wr <- sqrt(abs(term$gam$weights))
  w_resid <- term$gam$residuals * wr
  # Check that variances are actually available
  if (term$gam$Vp[1, 1] < 0) {
    stop("No variance estimates available")
  }
  list(
    covariance = covariance,
    w_resid = w_resid
  )
}
