.get_fit_plot_data <- function(
    pred_matrix_and_aux,
    term,
    gam) {
  mgcv_term <- gam$smooth[[term$term_idx]]
  first <- mgcv_term$first.para
  last <- mgcv_term$last.para
  attr(mgcv_term, "coefficients") <- gam$coefficients[first:last] # Relevant coeffs for i-th smooth

  pred_matrix <- pred_matrix_and_aux$X
  aux <- pred_matrix_and_aux$aux
  aux$smooth <- mgcv_term
  p <- gam$coefficients[first:last] ## relevant coefficients
  offset <- attr(pred_matrix, "offset") ## any term specific offset
  ## get fitted values ....
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
    pred_and_aux,
    se_with_mean) {
  first <- mgcv_term$first.para
  last <- mgcv_term$last.para
  attr(mgcv_term, "coefficients") <- gam$coefficients[first:last]

  if (se_with_mean && attr(mgcv_term, "nCons") > 0) {
    if (length(gam$cmX) < ncol(gam$Vp)) {
      gam$cmX <- c(gam$cmX, rep(0, ncol(gam$Vp) - length(gam$cmX))) # TODO does the mutation persist?
    }
    X1 <- matrix(gam$cmX, nrow(pred_and_aux$X), ncol(gam$Vp), byrow = TRUE)
    meanL1 <- mgcv_term$meanL1
    if (!is.null(meanL1)) {
      X1 <- X1 / meanL1
    }
    X1[, first:last] <- pred_and_aux$X
    se_fit <- sqrt(pmax(0, rowSums((X1 %*% gam$Vp) * X1)))
  } else {
    se_fit <- sqrt(pmax(0, rowSums((pred_and_aux$X %*% gam$Vp[first:last, first:last, drop = FALSE]) * pred_and_aux$X)))
  }
  se_fit
}
