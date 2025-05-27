.get_fit_and_errors_plot_data <- function(
    term,
    partial_resids,
    se,
    n,
    n2,
    ylim,
    xlim,
    too_far,
    se_with_mean,
    fit_smooth,
    w_resid,
    res_den,
    nsim,
    ...) {
  gam_viz <- term$gam_viz
  mgcv_term <- term$gam_viz$smooth[[term$term_idx]]
  first <- mgcv_term$first.para
  last <- mgcv_term$last.para
  attr(mgcv_term, "coefficients") <- gam_viz$coefficients[first:last] # Relevant coeffs for i-th smooth

  predict_matrix_and_aux <- .get_plot_predict_matrix_and_aux(
    mgcv_term = mgcv_term,
    data = gam_viz$model,
    n = n,
    n2 = n2,
    ylim = ylim,
    xlim = xlim,
    too_far = too_far,
    ...
  )
  if (is.null(predict_matrix_and_aux)) {
    stop("Cannot plot this term type.")
  } else {
    predict_matrix_and_aux$smooth <- mgcv_term
    if (is.null(predict_matrix_and_aux$fit)) {
      p <- gam_viz$coefficients[first:last] ## relevant coefficients
      offset <- attr(predict_matrix_and_aux$X, "offset") ## any term specific offset
      ## get fitted values ....
      if (is.null(offset)) {
        predict_matrix_and_aux$fit <- predict_matrix_and_aux$X %*% p
      } else {
        predict_matrix_and_aux$fit <- predict_matrix_and_aux$X %*% p + offset
      }
      if (!is.null(predict_matrix_and_aux$exclude)) {
        predict_matrix_and_aux$fit[predict_matrix_and_aux$exclude] <- NA
      }
      if (se) {
        ## get standard errors for fit

        ## test whether mean variability to be added to variability (only for centred terms)
        if (se_with_mean && attr(mgcv_term, "nCons") > 0) {
          if (length(gam_viz$cmX) < ncol(gam_viz$Vp)) {
            gam_viz$cmX <- c(gam_viz$cmX, rep(0, ncol(gam_viz$Vp) - length(gam_viz$cmX)))
          }
          X1 <- matrix(gam_viz$cmX, nrow(predict_matrix_and_aux$X), ncol(gam_viz$Vp), byrow = TRUE)
          meanL1 <- mgcv_term$meanL1
          if (!is.null(meanL1)) {
            X1 <- X1 / meanL1
          }
          X1[, first:last] <- predict_matrix_and_aux$X
          if (nsim > 0) {
            predict_matrix_and_aux$simF <- drop(predict_matrix_and_aux$fit) +
              X1 %*% t(rmvn(nsim, numeric(ncol(gam_viz$Vp)), gam_viz$Vp))
          }
          se.fit <- sqrt(pmax(0, rowSums((X1 %*% gam_viz$Vp) * X1)))
        } else {
          ## se in centred (or anyway unconstrained) space only
          if (nsim > 0) {
            predict_matrix_and_aux$simF <- drop(predict_matrix_and_aux$fit) +
              predict_matrix_and_aux$X %*%
              t(rmvn(
                nsim,
                numeric(length(p)),
                gam_viz$Vp[first:last, first:last, drop = FALSE]
              ))
          }
          se.fit <- sqrt(pmax(
            0,
            rowSums((predict_matrix_and_aux$X %*% gam_viz$Vp[first:last, first:last, drop = FALSE]) * predict_matrix_and_aux$X)
          ))
        }
        if (!is.null(predict_matrix_and_aux$exclude)) {
          predict_matrix_and_aux$se.fit[predict_matrix_and_aux$exclude] <- NA
        }
      } ## standard errors for fit completed
      if (partial_resids || (res_den != "none")) {
        predict_matrix_and_aux$p.resid <- fit_smooth + w_resid
      }
      if (se) {
        predict_matrix_and_aux$se <- se.fit
      }
      predict_matrix_and_aux$X <- NULL
    } else {
      ## P$fit created directly
      if (partial_resids || (res_den != "none")) {
        predict_matrix_and_aux$p.resid <- fit_smooth + w_resid
      }
    }
    predict_matrix_and_aux$plot.me <- TRUE
  }
  return(predict_matrix_and_aux)
}
