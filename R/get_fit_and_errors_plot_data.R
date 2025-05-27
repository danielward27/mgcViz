.get_fit_and_errors_plot_data <- function(
    term,
    gam_viz,
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
  first <- term$first.para
  last <- term$last.para
  attr(term, "coefficients") <- gam_viz$coefficients[first:last] # Relevant coeffs for i-th smooth

  P <- .get_plot_prediction_matrix_and_aux(
    term = term,
    data = gam_viz$model,
    n = n,
    n2 = n2,
    ylim = ylim,
    xlim = xlim,
    too_far = too_far,
    ...
  )
  if (is.null(P)) {
    P <- list(plot.me = FALSE)
  } else {
    P$smooth <- term
    if (is.null(P$fit)) {
      p <- gam_viz$coefficients[first:last] ## relevant coefficients
      offset <- attr(P$X, "offset") ## any term specific offset
      ## get fitted values ....
      if (is.null(offset)) {
        P$fit <- P$X %*% p
      } else {
        P$fit <- P$X %*% p + offset
      }
      if (!is.null(P$exclude)) {
        P$fit[P$exclude] <- NA
      }
      if (se) {
        ## get standard errors for fit

        ## test whether mean variability to be added to variability (only for centred terms)
        if (se_with_mean && attr(term, "nCons") > 0) {
          if (length(gam_viz$cmX) < ncol(gam_viz$Vp)) {
            gam_viz$cmX <- c(gam_viz$cmX, rep(0, ncol(gam_viz$Vp) - length(gam_viz$cmX)))
          }
          X1 <- matrix(gam_viz$cmX, nrow(P$X), ncol(gam_viz$Vp), byrow = TRUE)
          meanL1 <- term$meanL1
          if (!is.null(meanL1)) {
            X1 <- X1 / meanL1
          }
          X1[, first:last] <- P$X
          if (nsim > 0) {
            P$simF <- drop(P$fit) +
              X1 %*% t(rmvn(nsim, numeric(ncol(gam_viz$Vp)), gam_viz$Vp))
          }
          se.fit <- sqrt(pmax(0, rowSums((X1 %*% gam_viz$Vp) * X1)))
        } else {
          ## se in centred (or anyway unconstrained) space only
          if (nsim > 0) {
            P$simF <- drop(P$fit) +
              P$X %*%
              t(rmvn(
                nsim,
                numeric(length(p)),
                gam_viz$Vp[first:last, first:last, drop = FALSE]
              ))
          }
          se.fit <- sqrt(pmax(
            0,
            rowSums((P$X %*% gam_viz$Vp[first:last, first:last, drop = FALSE]) * P$X)
          ))
        }
        if (!is.null(P$exclude)) {
          P$se.fit[P$exclude] <- NA
        }
      } ## standard errors for fit completed
      if (partial_resids || (res_den != "none")) {
        P$p.resid <- fit_smooth + w_resid
      }
      if (se) {
        P$se <- se.fit
      }
      P$X <- NULL
    } else {
      ## P$fit created directly
      if (partial_resids || (res_den != "none")) {
        P$p.resid <- fit_smooth + w_resid
      }
    }
    P$plot.me <- TRUE
  }
  return(P)
}
