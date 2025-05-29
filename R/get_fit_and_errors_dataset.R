.get_covariance_and_errors <- function(
    term,
    gam,
    unconditional) { # TODO still doing two relatively unrelated actions...
  # Use Bayesian cov matrix including smoothing parameter uncertainty?
  covariance <- NULL
  if (unconditional) {
    if (is.null(gam$Vc)) {
      warning("Smoothness uncertainty corrected covariance not available")
    } else {
      covariance <- gam$Vc
    }
  }

  # produce working residuals if info available
  if (is.null(gam$residuals) || is.null(gam$weights)) {
    stop("Cannot compute working residuals.")
  }
  wr <- sqrt(abs(gam$weights))
  w_resid <- gam$residuals * wr
  # Check that variances are actually available
  if (gam$Vp[1, 1] < 0) {
    stop("No variance estimates available")
  }
  list(
    covariance = covariance,
    w_resid = w_resid
  )
}
