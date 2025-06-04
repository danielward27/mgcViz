.get_fit_plot_data <- function(
    pred_matrix,
    term,
    gam) {
  mgcv_term <- gam$smooth[[term$term_idx]]
  first <- mgcv_term$first.para
  last <- mgcv_term$last.para

  p <- gam$coefficients[first:last]
  offset <- attr(pred_matrix, "offset")

  if (is.null(offset)) {
    fit <- pred_matrix %*% p
  } else {
    fit <- pred_matrix %*% p + offset
  }
  fit
}

.get_errors_plot_data <- function(
    gam,
    mgcv_term,
    pred_matrix,
    se_with_mean) {
  first <- mgcv_term$first.para
  last <- mgcv_term$last.para
  if (se_with_mean && attr(mgcv_term, "nCons") > 0) {
    cmX <- gam$cmX
    if (length(cmX) < ncol(gam$Vp)) {
      cmX <- c(cmX, rep(0, ncol(gam$Vp) - length(cmX)))
    }
    X1 <- matrix(cmX, nrow(pred_matrix), ncol(gam$Vp), byrow = TRUE)
    meanL1 <- mgcv_term$meanL1
    if (!is.null(meanL1)) {
      X1 <- X1 / meanL1
    }
    X1[, first:last] <- pred_matrix
    se_fit <- sqrt(pmax(0, rowSums((X1 %*% gam$Vp) * X1)))
  } else {
    se_fit <- sqrt(
      pmax(
        0,
        rowSums((pred_matrix %*% gam$Vp[first:last, first:last, drop = FALSE]) * pred_matrix
        )))
  }
  se_fit
}
