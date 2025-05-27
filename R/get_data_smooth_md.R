#' @description Get data for plotting a 2D slice of a higher-dimensional smooth effects.
#' @export
#'
get_data.mgcv.smooth.MD <- function(
    term,
    fix,
    n = 40,
    xlim = NULL,
    ylim = NULL,
    maxpo = 1e4,
    too_far = c(0.1, NA),
    trans = identity,
    seWithMean = FALSE,
    unconditional = FALSE,
    ...) {
  if (length(too_far) == 1) {
    too_far <- c(too_far, NA)
  }

  P <- prepareP(
    term = term,
    unconditional = unconditional,
    residuals = TRUE,
    resDen = "none",
    se = TRUE,
    se_mult = 1,
    n = NULL,
    n2 = n,
    ylim = ylim,
    xlim = xlim,
    too_far = too_far,
    seWithMean = seWithMean,
    fix = fix
  )

  out <- .get_data_shared_2d(term = P$smooth, P = P, trans = trans, maxpo = maxpo)
  return(out)
}



# Internal function for preparing plot of two dimensional smooths
.prepare_plot_smooth_md <- function(
    term,
    fix,
    data = NULL,
    se_mult = 2,
    n2 = 40,
    ylim = NULL,
    xlim = NULL,
    too_far = 0.1,
    ...) {
  out <- NULL
  if (term$plot.me) {
    ov <- names(fix)
    iv <- term$term[!(term$term %in% ov)]
    xterm <- iv[1]
    yterm <- iv[2]
    raw <- data.frame(
      x = as.numeric(data[xterm][[1]]),
      y = as.numeric(data[yterm][[1]])
    )
    n2 <- max(10, n2)
    if (is.null(xlim)) {
      xm <- seq(min(raw$x), max(raw$x), length = n2)
    } else {
      xm <- seq(xlim[1], xlim[2], length = n2)
    }
    if (is.null(ylim)) {
      ym <- seq(min(raw$y), max(raw$y), length = n2)
    } else {
      ym <- seq(ylim[1], ylim[2], length = n2)
    }
    xx <- rep(xm, n2)
    yy <- rep(ym, each = n2)

    # Mark cells on X-Y grid are too far from any observation
    if (too_far[1] > 0) {
      exclude <- exclude.too.far(xx, yy, raw$x, raw$y, dist = too_far[1])
    } else {
      exclude <- rep(FALSE, n2 * n2)
    }

    # Mark covariate vectors (and corresponding residuals) that are too
    # far from X-Y plane (the slice of interest)
    if (is.na(too_far[2]) || too_far[2] > 0) {
      tmp <- sapply(ov, function(.nm) as.numeric(data[.nm][[1]]))
      tmp <- sqrt(mahalanobis(tmp, fix, diag(diag(cov(tmp)), ncol(tmp)))) # Euclidean distance
      exclude2 <- tmp >
        if (is.na(too_far[2])) {
          quantile(tmp, 0.1)
        } else {
          too_far[2]
        }
    } else {
      exclude2 <- FALSE
    }
    if (term$by != "NA") {
      # deal with any by variables
      by <- rep(1, n2^2)
      dat <- data.frame(x = xx, y = yy, by = by)
      colnames(dat) <- c(xterm, yterm, term$by)
    } else {
      dat <- data.frame(x = xx, y = yy)
      colnames(dat) <- c(xterm, yterm)
    } ## prediction data.frame complete

    for (ii in ov) {
      dat[[ii]] <- rep(fix[ii], n2^2)
    }

    X <- PredictMat(term, dat) ## prediction matrix for this term

    if (is.null(ylim)) {
      ylim <- range(ym)
    }
    if (is.null(xlim)) {
      xlim <- range(xm)
    }
    out <- list(
      X = X,
      x = xm,
      y = ym,
      scale = FALSE,
      se = TRUE,
      raw = raw,
      se_mult = se_mult,
      ylim = ylim,
      xlim = xlim,
      exclude = exclude,
      exclude2 = exclude2
    )
  }
  return(out)
}
