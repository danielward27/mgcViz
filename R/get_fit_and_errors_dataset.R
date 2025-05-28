.get_covariance_and_errors <- function(
    term,
    unconditional) { # TODO still doing two relatively unrelated actions...
  # Use Bayesian cov matrix including smoothing parameter uncertainty?
  covariance <- NULL
  if (unconditional) {
    if (is.null(term$gam_viz$Vc)) {
      warning("Smoothness uncertainty corrected covariance not available")
    } else {
      covariance <- term$gam_viz$Vc
    }
  }

  # TODO more descriptive name for residuals, looks like below we get Pearson residuals,
  # not absolute.

  # produce working residuals if info available
  if (is.null(term$gam_viz$residuals) || is.null(term$gam_viz$weights)) {
    stop("Cannot compute working residuals.")
  }
  wr <- sqrt(abs(term$gam_viz$weights))
  w_resid <- term$gam_viz$residuals * wr
  # Check that variances are actually available
  if (term$gam_viz$Vp[1, 1] < 0) {
    stop("No variance estimates available")
  }
  list(
    covariance = covariance,
    w_resid = w_resid
  )
}
