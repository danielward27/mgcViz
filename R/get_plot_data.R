#' @param sm XXX
#' @param x  XXX
#' @param partial_resids  XXX
#' @param se  XXX
#' @param n  XXX
#' @param n2  XXX
#' @param ylim  XXX
#' @param xlim  XXX
#' @param too_far  XXX
#' @param se1_mult  XXX
#' @param se2_mult  XXX
#' @param seWithMean  XXX
#' @param fitSmooth  XXX
#' @param w_resid  XXX
#' @param resDen  XXX
#' @param ...  XXX
#' @noRd
#' @examples
#' library(mgcViz)
#' n <- 1e3
#' x1 <- rnorm(n)
#' x2 <- rnorm(n)
#' dat <- data.frame(
#'   "x1" = x1, "x2" = x2,
#'   "y" = sin(x1) + 0.5 * x2^2 + pmax(x2, 0.2) * rnorm(n)
#' )
#' b <- bam(y ~ s(x1) + s(x2), data = dat, method = "fREML", discrete = TRUE)
#' v <- getViz(b)(1)
#' class(v)
#' o <- v
#' o$smooth <- o$gam_viz$smooth[[o$ism]]
#' fv_terms <- o$store$termsFit[, o$store$np + o$ism]
#' init <- mgcViz:::.initializeXXX(o, unconditional = FALSE, residuals = FALSE, resDen = "cond", se = TRUE, fv_terms)
#' o <- init$o
#' w_resid <- init$w_resid
#' partial_resids <- init$partial_resids
#' se2_mult <- init$se2_mult
#' se1_mult <- init$se1_mult
#' se <- init$se
#' fv_terms <- init$fv_terms
#' order <- init$order
#' sm <- o$smooth
#' too_far <- 0.1
#' seWithMean <- FALSE
#' resDen <- "none"
#' @noRd
.get_plot_data <- function(
    term,
    gam_viz,
    partial_resids,
    se,
    n,
    n2,
    ylim,
    xlim,
    too_far,
    se1_mult,
    se2_mult,
    seWithMean,
    fitSmooth,
    w_resid,
    resDen,
    nsim,
    ...) {
  first <- term$first.para
  last <- term$last.para
  edf <- sum(gam_viz$edf[first:last]) ## Effective DoF for this term
  attr(term, "coefficients") <- gam_viz$coefficients[first:last] # Relevant coeffs for i-th smooth
  P <- .prepare(
    term = term,
    data = gam_viz$model,
    se1_mult = se1_mult,
    se2_mult = se2_mult,
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
      if (se && P$se) {
        ## get standard errors for fit
        ## test whether mean variability to be added to variability (only for centred terms)
        if (seWithMean && attr(term, "nCons") > 0) {
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
          se.fit <- sqrt(pmax(0, rowSums((X1 %*%gam_viz$Vp) * X1)))
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
      if (partial_resids || (resDen != "none")) {
        P$p.resid <- fitSmooth + w_resid
      }
      if (se && P$se) {
        P$se <- se.fit * P$se_mult
      } # Note multiplier
      P$X <- NULL
    } else {
      ## P$fit created directly
      if (partial_resids || (resDen != "none")) {
        P$p.resid <- fitSmooth + w_resid
      }
    }
    P$plot.me <- TRUE
  }
  return(P)
}
