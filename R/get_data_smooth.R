#' @description Get data for plotting one dimensional smooth effects.
#' @export
get_data.mgcv.smooth.1D <- function(
    term,
    n = 100,
    xlim = NULL,
    maxpo = 1e4,
    trans = identity,
    unconditional = FALSE,
    se_with_mean = FALSE,
    nsim = 0,
    ...) {
  P <- prepareP(
    term = term,
    unconditional = unconditional,
    residuals = TRUE,
    res_den = "none",
    n = n,
    n2 = NULL,
    ylim = NULL,
    xlim = xlim,
    too_far = NULL,
    se_with_mean = se_with_mean,
    nsim = nsim
  )
  .dat <- list()
  if (!is.null(P$raw)) {
    # Construct data.frame of partial residuals
    res <- data.frame("x" = as.vector(P$raw))
    if (!is.null(P$p_resid) & length(P$p_resid)) {
      res$y <- trans(P$p_resid)
    }

    # Exclude residuals falling outside boundaries
    .dat$res <- res[res$x >= P$xlim[1] & res$x <= P$xlim[2], , drop = FALSE]

    # Sample if too many points (> maxpo)
    nres <- nrow(.dat$res)
    .dat$res$sub <- if (nres > maxpo) {
      sample(c(rep(T, maxpo), rep(F, nres - maxpo)))
    } else {
      rep(T, nres)
    }
  }

  .dat$fit <- data.frame(
    x = P$x, # x values
    y = P$fit, # fitted values
    ty = trans(P$fit), # fitted values after trans
    se = P$se
  ) # standard error

  if (!is.null(P$sim_f)) {
    nsim <- ncol(P$sim_f)
    .dat$sim <- data.frame(
      "x" = rep(P$x, nsim),
      "ty" = trans(as.vector(P$sim_f)),
      "id" = as.factor(rep(1:nsim, each = nrow(P$sim_f)))
    )
  }

  .dat$misc <- list("trans" = trans)
  return(.dat)
}



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
    se_with_mean = FALSE,
    unconditional = FALSE,
    ...) {
  # 1) Prepare data
  P <- prepareP(
    term = term,
    unconditional = unconditional,
    residuals = TRUE,
    res_den = "none",
    n = NULL,
    n2 = n,
    ylim = ylim,
    xlim = xlim,
    too_far = too_far,
    se_with_mean = se_with_mean
  )

  # 2) Produce output object
  out <- .get_data_shared_2d(P = P, trans = trans, maxpo = maxpo)
  return(out)
}

# TODO do we need both generics for get_data and

# Used by e.g MD and SOS
#' @noRd
.get_data_shared_2d <- function(P, trans, maxpo, flip = FALSE) {
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
    se_with_mean = FALSE,
    unconditional = FALSE,
    ...) {
  if (length(too_far) == 1) {
    too_far <- c(too_far, NA)
  }

  P <- prepareP(
    term = term,
    unconditional = unconditional,
    residuals = TRUE,
    res_den = "none",
    n = NULL,
    n2 = n,
    ylim = ylim,
    xlim = xlim,
    too_far = too_far,
    se_with_mean = se_with_mean,
    fix = fix
  )

  out <- .get_data_shared_2d(P = P, trans = trans, maxpo = maxpo)
  return(out)
}



#' @description Default plot preparation method for smooth objects `x'
#' inheriting from "mgcv.smooth".
#' @export
.get_plot_prediction_matrix_and_aux.mgcv.smooth <- function(
    mgcv_term,
    data = NULL,
    n = 100,
    n2 = 40,
    ylim = NULL,
    xlim = NULL,
    too_far = 0.1,
    ...) {
  if (mgcv_term$dim == 1) {
    out <- .get_plot_prediction_matrix_and_aux_plot_smooth_1d(
      mgcv_term = mgcv_term,
      data = data,
      n = n,
      xlim = xlim,
      ...
    )
  }

  if (mgcv_term$dim == 2) {
    out <- .get_plot_prediction_matrix_and_aux_plot_smooth_2d(
      mgcv_term = mgcv_term,
      data = data,
      n2 = n2,
      ylim = ylim,
      xlim = xlim,
      too_far = too_far,
      ...
    )
  }

  if (mgcv_term$dim > 2) {
    out <- .get_plot_prediction_matrix_and_aux_plot_smooth_md(
      mgcv_term = mgcv_term,
      data = data,
      n2 = n2,
      ylim = ylim,
      xlim = xlim,
      too_far = too_far,
      ...
    )
  }

  return(out)
}


# Internal function for preparing plot of one dimensional smooths
.get_plot_prediction_matrix_and_aux_plot_smooth_1d <- function(
    mgcv_term,
    data,
    n = 100,
    xlim = NULL,
    ...) {
  raw <- as.vector(data[mgcv_term$term][[1]])
  if (is.null(xlim)) {
    # Generate x sequence for prediction
    x_seq <- seq(min(raw), max(raw), length = n)
  } else {
    x_seq <- seq(xlim[1], xlim[2], length = n)
  }
  if (mgcv_term$by != "NA") {
    # Deal with any by variables
    by <- rep(1, n)
    dat <- data.frame(x = x_seq, by = by)
    names(dat) <- c(mgcv_term$term, mgcv_term$by)
  } else {
    dat <- data.frame(x = x_seq)
    names(dat) <- mgcv_term$term
  } # Finished preparing prediction data.frame
  X <- PredictMat(mgcv_term, dat) # prediction matrix for this term

  if (is.null(xlim)) {
    xlim <- range(x_seq)
  }
  out <- list(
    X = X,
    x = x_seq,
    raw = raw
  )
  return(out)
}

# Internal function for preparing plot of two dimensional smooths
.get_plot_prediction_matrix_and_aux_plot_smooth_2d <- function(
    mgcv_term,
    data = NULL,
    n2 = 40,
    ylim = NULL, # TODO name y is confusing!
    xlim = NULL,
    too_far = 0.1,
    ...) {
  xterm <- mgcv_term$term[1]
  yterm <- mgcv_term$term[2]
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
  if (mgcv_term$by != "NA") {
    # deal with any by variables
    by <- rep(1, n2^2)
    dat <- data.frame(x = x_rep, y = y_rep, by = by)
    colnames(dat) <- c(xterm, yterm, mgcv_term$by)
  } else {
    dat <- data.frame(x = x_rep, y = y_rep)
    colnames(dat) <- c(xterm, yterm)
  } ## prediction data.frame complete
  X <- PredictMat(mgcv_term, dat) ## prediction matrix for this term

  out <- list(
    X = X,
    x = x_seq,
    y = y_seq,
    raw = raw,
    exclude = exclude
  )

  return(out)
}


# Internal function for preparing plot of two dimensional smooths
.get_plot_prediction_matrix_and_aux_plot_smooth_md <- function(
    mgcv_term,
    fix,
    data = NULL,
    n2 = 40,
    ylim = NULL,
    xlim = NULL,
    too_far = 0.1,
    ...) {
  ov <- names(fix)
  iv <- mgcv_term$term[!(mgcv_term$term %in% ov)]
  xterm <- iv[1]
  yterm <- iv[2]
  raw <- data.frame(
    x = as.numeric(data[xterm][[1]]),
    y = as.numeric(data[yterm][[1]])
  )
  n2 <- max(10, n2)
  if (is.null(xlim)) {
    x_seq <- seq(min(raw$x), max(raw$x), length = n2)
  } else {
    x_seq <- seq(xlim[1], xlim[2], length = n2)
  }
  if (is.null(ylim)) {
    y_seq <- seq(min(raw$y), max(raw$y), length = n2)
  } else {
    y_seq <- seq(ylim[1], ylim[2], length = n2)
  }
  xx <- rep(x_seq, n2)
  yy <- rep(y_seq, each = n2)

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
  if (mgcv_term$by != "NA") {
    # deal with any by variables
    by <- rep(1, n2^2)
    dat <- data.frame(x = xx, y = yy, by = by)
    colnames(dat) <- c(xterm, yterm, mgcv_term$by)
  } else {
    dat <- data.frame(x = xx, y = yy)
    colnames(dat) <- c(xterm, yterm)
  } ## prediction data.frame complete

  for (ii in ov) {
    dat[[ii]] <- rep(fix[ii], n2^2)
  }

  X <- PredictMat(mgcv_term, dat) ## prediction matrix for this term

  if (is.null(ylim)) {
    ylim <- range(y_seq)
  }
  if (is.null(xlim)) {
    xlim <- range(x_seq)
  }
  out <- list(
    X = X,
    x = x_seq,
    y = y_seq,
    raw = raw,
    exclude = exclude,
    exclude2 = exclude2
  )
  return(out)
}
