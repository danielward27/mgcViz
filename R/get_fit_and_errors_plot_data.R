.get_fit_and_errors_plot_data <- function(
    pred_matrix_and_aux,
    term,
    gam,
    compute_partial_resids,
    compute_se,
    se_with_mean,
    term_fit,
    w_resid,
    nsim) {
  gam <- gam
  mgcv_term <- gam$smooth[[term$term_idx]]
  first <- mgcv_term$first.para
  last <- mgcv_term$last.para
  attr(mgcv_term, "coefficients") <- gam$coefficients[first:last] # Relevant coeffs for i-th smooth

  pred_matrix <- pred_matrix_and_aux$X
  aux <- pred_matrix_and_aux$aux

  if (is.null(pred_matrix_and_aux)) {
    stop("Cannot plot this term type.")
  } else {
    aux$smooth <- mgcv_term
    p <- gam$coefficients[first:last] ## relevant coefficients
    offset <- attr(pred_matrix, "offset") ## any term specific offset
    ## get fitted values ....
    if (is.null(offset)) {
      fit <- pred_matrix %*% p
    } else {
      fit <- pred_matrix %*% p + offset
    }
    if (!is.null(pred_matrix_and_aux$exclude)) {
      fit[pred_matrix_and_aux$exclude] <- NA
    }
    if (compute_se) {
      ## get standard errors for fit

      ## test whether mean variability to be added to variability (only for centred terms)
      if (se_with_mean && attr(mgcv_term, "nCons") > 0) {
        if (length(gam$cmX) < ncol(gam$Vp)) {
          gam$cmX <- c(gam$cmX, rep(0, ncol(gam$Vp) - length(gam$cmX)))
        }
        X1 <- matrix(gam$cmX, nrow(pred_matrix), ncol(gam$Vp), byrow = TRUE)
        meanL1 <- mgcv_term$meanL1
        if (!is.null(meanL1)) {
          X1 <- X1 / meanL1
        }
        X1[, first:last] <- pred_matrix
        if (nsim > 0) {
          pred_matrix_and_aux$simF <- drop(fit) +
            X1 %*% t(rmvn(nsim, numeric(ncol(gam$Vp)), gam$Vp))
        }
        se.fit <- sqrt(pmax(0, rowSums((X1 %*% gam$Vp) * X1)))
      } else {
        ## se in centred (or anyway unconstrained) space only
        if (nsim > 0) {
          pred_matrix_and_aux$simF <- drop(fit) +
            pred_matrix %*%
            t(rmvn(
              nsim,
              numeric(length(p)),
              gam$Vp[first:last, first:last, drop = FALSE]
            ))
        }
        se.fit <- sqrt(pmax(
          0,
          rowSums((pred_matrix %*% gam$Vp[first:last, first:last, drop = FALSE]) * pred_matrix)
        ))
      }
      if (!is.null(pred_matrix_and_aux$exclude)) {
        se.fit[pred_matrix_and_aux$exclude] <- NA
      }
    } ## standard errors for fit completed
    if (compute_partial_resids) {
      partial_resids <- term_fit + w_resid
    }
  }
  list(
    fit = fit,
    partial_resids = partial_resids,
    se.fit = se.fit
  )
}
