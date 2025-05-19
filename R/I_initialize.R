#' .initializeXXX
#'
#' @param o, ...
#' @param unconditional, ...
#' @param residuals, ...
#' @param resDen, ...
#' @param se, ...
#' @param fv.terms, ...
#' @return a list
#' @noRd
.initialize <- function(term, unconditional, residuals, resDen, se, se.mult) {
  V <- fv.terms <- NULL

  # Use Bayesian cov matrix including smoothing parameter uncertainty?
  if (unconditional) {
    if (is.null(term$gObj$Vc)) {
      warning("Smoothness uncertainty corrected covariance not available")
    } else {
      V <- term$gObj$Vc
    }
  }

  w.resid <- NULL
  if (length(residuals) > 1) {
    # residuals supplied
    if (length(residuals) == length(term$gObj$residuals)) {
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
      if (is.null(term$gObj$residuals) || is.null(term$gObj$weights)) {
        partial.resids <- FALSE
      } else {
        wr <- sqrt(abs(term$gObj$weights))
        w.resid <- term$gObj$residuals * wr
      }
    }

    fv.terms <- term$gObj$store$termsFit[, term$gObj$store$np + term$ism]
    if (is.null(fv.terms)) {
      fv.terms <- predict(term$gObj, type = "terms")
    }
  }

  if (se) {
    # Sort out CI widths for 1D and 2D smooths
    if (se.mult < 0) {
      se.mult <- 0
    }
    # Check that variances are actually available
    if (term$gObj$Vp[1, 1] < 0) {
      se <- FALSE
      warning("No variance estimates available")
    }
  }

  ## Array giving the order of each parametric term
  order <- if (is.list(term$gObj$pterms)) {
    unlist(lapply(term$gObj$pterms, attr, "order"))
  } else {
    attr(term$gObj$pterms, "order")
  }

  return(list(
    V = V,
    order = order,
    w.resid = w.resid,
    se.mult = se.mult,
    se = se,
    fv.terms = fv.terms,
    partial.resids = partial.resids
  ))
}
