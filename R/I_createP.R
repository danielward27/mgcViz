#' @param sm XXX
#' @param x  XXX
#' @param partial.resids  XXX
#' @param se  XXX
#' @param n  XXX
#' @param n2  XXX
#' @param ylim  XXX
#' @param xlim  XXX
#' @param too.far  XXX
#' @param se1.mult  XXX
#' @param se2.mult  XXX
#' @param seWithMean  XXX
#' @param fitSmooth  XXX
#' @param w.resid  XXX
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
#' o$smooth <- o$gObj$smooth[[o$ism]]
#' fv.terms <- o$store$termsFit[, o$store$np + o$ism]
#' init <- mgcViz:::.initializeXXX(o, unconditional = FALSE, residuals = FALSE, resDen = "cond", se = TRUE, fv.terms)
#' o <- init$o
#' w.resid <- init$w.resid
#' partial.resids <- init$partial.resids
#' se2.mult <- init$se2.mult
#' se1.mult <- init$se1.mult
#' se <- init$se
#' fv.terms <- init$fv.terms
#' order <- init$order
#' sm <- o$smooth
#' too.far <- 0.1
#' seWithMean <- FALSE
#' resDen <- "none"
#' @noRd
.createP <- function(
    sm,
    gam, # I think this is an mgcv GAM
    partial.resids,
    se,
    n,
    n2,
    ylim,
    xlim,
    too.far,
    se1.mult,
    se2.mult,
    seWithMean,
    fitSmooth,
    w.resid,
    resDen,
    nsim,
    ...) {
  first <- sm$first.para
  last <- sm$last.para
  edf <- sum(gam$edf[first:last]) ## Effective DoF for this term
  attr(sm, "coefficients") <- gam$coefficients[first:last] # Relevant coeffs for i-th smooth
  P <- .prepare(
    x = sm,
    data = gam$model,
    se1.mult = se1.mult,
    se2.mult = se2.mult,
    n = n,
    n2 = n2,
    ylim = ylim,
    xlim = xlim,
    too.far = too.far,
    ...
  )
  if (is.null(P)) {
    P <- list(plot.me = FALSE)
  } else {
    P$smooth <- sm
    if (is.null(P$fit)) {
      p <- gam$coefficients[first:last] ## relevant coefficients
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
        if (seWithMean && attr(sm, "nCons") > 0) {
          if (length(gam$cmX) < ncol(gam$Vp)) {
            gam$cmX <- c(gam$cmX, rep(0, ncol(gam$Vp) - length(gam$cmX)))
          }
          X1 <- matrix(gam$cmX, nrow(P$X), ncol(gam$Vp), byrow = TRUE)
          meanL1 <- sm$meanL1
          if (!is.null(meanL1)) {
            X1 <- X1 / meanL1
          }
          X1[, first:last] <- P$X
          if (nsim > 0) {
            P$simF <- drop(P$fit) +
              X1 %*% t(rmvn(nsim, numeric(ncol(gam$Vp)), gam$Vp))
          }
          se.fit <- sqrt(pmax(0, rowSums((X1 %*%gam$Vp) * X1)))
        } else {
          ## se in centred (or anyway unconstrained) space only
          if (nsim > 0) {
            P$simF <- drop(P$fit) +
              P$X %*%
              t(rmvn(
                nsim,
                numeric(length(p)),
                gam$Vp[first:last, first:last, drop = FALSE]
              ))
          }
          se.fit <- sqrt(pmax(
            0,
            rowSums((P$X %*% gam$Vp[first:last, first:last, drop = FALSE]) * P$X)
          ))
        }
        if (!is.null(P$exclude)) {
          P$se.fit[P$exclude] <- NA
        }
      } ## standard errors for fit completed
      if (partial.resids || (resDen != "none")) {
        P$p.resid <- fitSmooth + w.resid
      }
      if (se && P$se) {
        P$se <- se.fit * P$se.mult
      } # Note multiplier
      P$X <- NULL
    } else {
      ## P$fit created directly
      if (partial.resids || (resDen != "none")) {
        P$p.resid <- fitSmooth + w.resid
      }
    }
    P$plot.me <- TRUE
  }
  return(P)
}
