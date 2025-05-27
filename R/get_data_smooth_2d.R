#' @description Get data for plotting two dimensional smooth effects.
#' @export
get_data.mgcv.smooth.2D <- function(
    term,
    n = 40,
    xlim = NULL,
    ylim = NULL,
    maxpo = 1e4,
    too_far = 0.1,
    trans = identity,
    seWithMean = FALSE,
    unconditional = FALSE,
    ...) {
  # 1) Prepare data
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
    seWithMean = seWithMean
  )

  # 2) Produce output object
  out <- .get_data_shared_2d(term = P$smooth, P = P, trans = trans, maxpo = maxpo)
  return(out)
}


# Used by e.g MD and SOS
#' @noRd
.get_data_shared_2d <- function(term, P, trans, maxpo, flip = FALSE) { # TODO is term even used?
  .dat <- list()
  # 1) Build dataset on fitted effect
  P$fit[P$exclude] <- NA
  .dat$fit <- data.frame(
    "z" = drop(P$fit),
    "tz" = drop(trans(P$fit)),
    "x" = rep(P$x, length(P$fit) / length(P$x)),
    "y" = rep(P$y, each = length(P$fit) / length(P$x)),
    "se" = P$se
  )

  # 2) Build dataset on residuals
  if (!is.null(P$raw)) {
    # Exclude points too far from current slice (relevant only when called by plot.mgcv.smooth.MD)
    if (!is.null(P$exclude2) && any(P$exclude2)) {
      P$raw <- P$raw[!P$exclude2, ]
    }

    # Exclude residuals falling outside boundaries
    .dat$res <- P$raw[
      P$raw$x >= P$xlim[1] &
        P$raw$x <= P$xlim[2] &
        P$raw$y >= P$ylim[1] &
        P$raw$y <= P$ylim[2], ,
      drop = FALSE
    ]

    # Sample if too many points (> maxpo)
    nres <- nrow(.dat$res)
    .dat$res$sub <- if (nres > maxpo) {
      sample(c(rep(T, maxpo), rep(F, nres - maxpo)))
    } else {
      rep(T, nres)
    }
  }

  .dat$misc <- list("trans" = trans)

  if (flip) {
    .dat2 <- .dat
    .dat2$fit$x <- .dat$fit$y
    .dat2$fit$y <- .dat$fit$x
    .dat2$res$x <- .dat$res$y
    .dat2$res$y <- .dat$res$x
    .dat <- .dat2
  }

  return(.dat)
}


# Internal function for preparing plot of two dimensional smooths
.prepare_plot_smooth_2d <- function(
    term,
    data = NULL,
    se_mult = 2,
    n2 = 40,
    ylim = NULL,
    xlim = NULL,
    too_far = 0.1,
    ...) {
  out <- NULL
  if (term$plot.me) {
    xterm <- term$term[1]
    yterm <- term$term[2]
    raw <- data.frame(
      x = as.numeric(data[xterm][[1]]),
      y = as.numeric(data[yterm][[1]])
    )
    n2 <- max(10, n2)


    if (is.null(xlim)) {
      xlim <- c(min(raw$x), max(raw$x))
    }
    if (is.null(ylim)) {
      ylim <- c(min(raw$y), max(raw$y))
    }

    x_seq <- seq(xlim[1], xlim[2], length = n2)
    y_seq <- seq(ylim[1], ylim[2], length = n2)
    x_rep <- rep(x_seq, n2)
    y_rep <- rep(y_seq, rep(n2, n2))
    if (too_far > 0) {
      exclude <- exclude.too.far(x_rep, y_rep, raw$x, raw$y, dist = too_far)
    } else {
      exclude <- rep(FALSE, n2 * n2)
    }
    if (term$by != "NA") {
      # deal with any by variables
      by <- rep(1, n2^2)
      dat <- data.frame(x = x_rep, y = y_rep, by = by)
      colnames(dat) <- c(xterm, yterm, term$by)
    } else {
      dat <- data.frame(x = x_rep, y = y_rep)
      colnames(dat) <- c(xterm, yterm)
    } ## prediction data.frame complete
    X <- PredictMat(term, dat) ## prediction matrix for this term

    out <- list(
      X = X,
      x = x_seq,
      y = y_seq,
      scale = FALSE,
      se = TRUE,
      raw = raw,
      se_mult = se_mult,
      ylim = ylim,
      xlim = xlim,
      exclude = exclude
    )
  }
  return(out)
}
