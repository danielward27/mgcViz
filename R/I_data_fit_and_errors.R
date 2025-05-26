
.data_fit_and_errors <- function(term, unconditional, residuals, resDen, se, se.mult) {
  V <- fv.terms <- NULL

  # Use Bayesian cov matrix including smoothing parameter uncertainty?
  if (unconditional) {
    if (is.null(term$gam_viz$Vc)) {
      warning("Smoothness uncertainty corrected covariance not available")
    } else {
      V <- term$gam_viz$Vc
    }
  }

  w.resid <- NULL
  if (length(residuals) > 1) {
    # residuals supplied
    if (length(residuals) == length(term$gam_viz$residuals)) {
      w.resid <- residuals
    } else {
      warning("residuals argument to plot.gamViz is wrong length: ignored")
    }
    partial.resids <- TRUE
  } else {
    partial.resids <- residuals
  } # use working residuals or none

  # Getting information needed for partial residuals
  if (partial.resids || (resDen != "none")) {
    if (is.null(w.resid)) {
      # produce working residuals if info available
      if (is.null(term$gam_viz$residuals) || is.null(term$gam_viz$weights)) {
        partial.resids <- FALSE
      } else {
        wr <- sqrt(abs(term$gam_viz$weights))
        w.resid <- term$gam_viz$residuals * wr
      }
    }

    fv.terms <- term$gam_viz$store$termsFit[, term$gam_viz$store$np + term$ism]
    if (is.null(fv.terms)) {
      fv.terms <- predict(term$gam_viz, type = "terms")
    }
  }

  if (se) {
    # Sort out CI widths for 1D and 2D smooths
    if (se.mult < 0) {
      se.mult <- 0
    }
    # Check that variances are actually available
    if (term$gam_viz$Vp[1, 1] < 0) {
      se <- FALSE
      warning("No variance estimates available")
    }
  }


  return(list(
    V = V,
    w.resid = w.resid,
    se.mult = se.mult,
    se = se,
    fv.terms = fv.terms,
    partial.resids = partial.resids
  ))
}
