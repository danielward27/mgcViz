.get_fit_and_errors_dataset <- function(term, unconditional, residuals, res_den) {
  V <- fv_terms <- NULL

  # Use Bayesian cov matrix including smoothing parameter uncertainty?
  if (unconditional) {
    if (is.null(term$gam_viz$Vc)) {
      warning("Smoothness uncertainty corrected covariance not available")
    } else {
      V <- term$gam_viz$Vc
    }
  }

  w_resid <- NULL
  if (length(residuals) > 1) {
    # residuals supplied
    if (length(residuals) == length(term$gam_viz$residuals)) {
      w_resid <- residuals
    } else {
      warning("residuals argument to plot.gamViz is wrong length: ignored")
    }
    partial_resids <- TRUE
  } else {
    partial_resids <- residuals
  } # use working residuals or none

  # Getting information needed for partial residuals
  if (partial_resids || (res_den != "none")) {
    if (is.null(w_resid)) {
      # produce working residuals if info available
      if (is.null(term$gam_viz$residuals) || is.null(term$gam_viz$weights)) {
        partial_resids <- FALSE
      } else {
        wr <- sqrt(abs(term$gam_viz$weights))
        w_resid <- term$gam_viz$residuals * wr
      }
    }

    fv_terms <- term$gam_viz$store$termsFit[, term$gam_viz$store$np + term$term_idx]
    if (is.null(fv_terms)) {
      fv_terms <- predict(term$gam_viz, type = "terms") # TODO inefficient?
    }
  }

  # Check that variances are actually available
  if (term$gam_viz$Vp[1, 1] < 0) {
    stop("No variance estimates available")
  }

  list(
    V = V,
    w_resid = w_resid,
    fv_terms = fv_terms,
    w_resid = w_resid,
    partial_resids = partial_resids,
    V = V
  )
}
