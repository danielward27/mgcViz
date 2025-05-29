#' @description Get data for plotting one dimensional smooth effects.
#' @export
get_data.mgcv.smooth.1D <- function(
    term,
    fitted_terms,
    n = 100,
    lims = NULL,
    maxpo = 1e4,
    trans = identity,
    unconditional = FALSE,
    se_with_mean = FALSE,
    nsim = 0,
    ...) {
  P <- prepareP(
    term = term,
    fitted_terms = fitted_terms,
    unconditional = unconditional,
    residuals = TRUE,
    n = n,
    n2 = NULL,
    lims = lims,
    too_far = NULL,
    se_with_mean = se_with_mean,
    nsim = nsim
  )
  .dat <- list() # TODO I don't like .dat
  if (!is.null(P$aux$raw)) {
    # Construct data.frame of partial residuals
    res <- data.frame("x" = as.vector(P$aux$raw))
    if (!is.null(P$partial_resids) && length(P$partial_resids)) {
      res$y <- trans(P$partial_resids)
    }
    # Exclude residuals falling outside boundaries
    .dat$res <- res
  }

  .dat$fit <- data.frame(
    x = P$aux$x, # x values
    y = P$fit, # fitted values
    ty = trans(P$fit), # fitted values after trans
    se = P$se
  ) # standard error

  if (!is.null(P$sim_f)) {
    nsim <- ncol(P$sim_f)
    .dat$sim <- data.frame(
      "x" = rep(P$aux$x, nsim),
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
    fitted_terms,
    n = 40,
    lims = NULL,
    too_far = 0.1,
    trans = identity,
    se_with_mean = FALSE,
    unconditional = FALSE,
    ...) {
  # 1) Prepare data
  P <- prepareP(
    term = term,
    fitted_terms = fitted_terms,
    unconditional = unconditional,
    residuals = TRUE,
    n = NULL,
    n2 = n,
    lims = lims,
    too_far = too_far,
    se_with_mean = se_with_mean
  )

  # 2) Produce output object
  out <- .get_data_shared_2d(P = P, trans = trans)
  return(out)
}

# TODO do we need both generics for get_data and

# Used by e.g MD and SOS
#' @noRd
.get_data_shared_2d <- function(P, trans, flip = FALSE) {
  .dat <- list()
  # 1) Build dataset on fitted effect
  P$fit[P$aux$exclude] <- NA
  .dat$fit <- data.frame(
    "z" = drop(P$fit),
    "tz" = drop(trans(P$fit)),
    "x" = rep(P$aux$x, length(P$fit) / length(P$aux$x)),
    "y" = rep(P$aux$y, each = length(P$fit) / length(P$aux$x)),
    "se" = P$se
  )

  .dat$res <- P$aux$raw
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
    fitted_terms,
    fix,
    n = 40,
    lims = NULL,
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
    fitted_terms = fitted_terms,
    unconditional = unconditional,
    residuals = TRUE,
    n = NULL,
    n2 = n,
    lims = lims,
    too_far = too_far,
    se_with_mean = se_with_mean,
    fix = fix
  )

  .get_data_shared_2d(P = P, trans = trans)
}



#' @description Default plot preparation method for smooth objects `x'
#' inheriting from "mgcv.smooth".
#' @export
.get_plot_predict_matrix_and_aux.mgcv.smooth <- function(
    mgcv_term,
    data = NULL,
    n = 100,
    n2 = 40,
    lims = NULL,
    too_far = 0.1,
    ...) {
  if (mgcv_term$dim == 1) {
    out <- .get_plot_predict_matrix_and_aux_plot_smooth_1d(
      mgcv_term = mgcv_term,
      data = data,
      n = n,
      lims = lims,
      ...
    )
  }

  if (mgcv_term$dim == 2) {
    out <- .get_plot_predict_matrix_and_aux_plot_smooth_2d(
      mgcv_term = mgcv_term,
      data = data,
      n2 = n2,
      lims = lims,
      too_far = too_far,
      ...
    )
  }

  if (mgcv_term$dim > 2) {
    out <- .get_plot_predict_matrix_and_aux_plot_smooth_md(
      mgcv_term = mgcv_term,
      data = data,
      n2 = n2,
      lims = lims,
      too_far = too_far,
      ...
    )
  }
  out
}


# Internal function for preparing plot of one dimensional smooths
.get_plot_predict_matrix_and_aux_plot_smooth_1d <- function(
    mgcv_term,
    data,
    n = 100,
    lims = NULL,
    ...) {
  raw <- as.vector(data[mgcv_term$term][[1]])
  if (is.null(lims)) {
    # Generate x sequence for prediction
    x_seq <- seq(min(raw), max(raw), length = n)
  } else {
    x_seq <- seq(lims[1], lims[2], length = n)
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

  out <- list(
    X = X,
    aux = list(x = x_seq, raw = raw)
  )
  return(out)
}

# Internal function for preparing plot of two dimensional smooths
.get_plot_predict_matrix_and_aux_plot_smooth_2d <- function(
    mgcv_term,
    data = NULL,
    n2 = 40,
    lims = NULL,
    too_far = 0.1,
    ...) {
  xterm <- mgcv_term$term[1]
  yterm <- mgcv_term$term[2]
  raw <- data.frame(
    x = as.numeric(data[xterm][[1]]),
    y = as.numeric(data[yterm][[1]])
  )
  if (is.null(lims)) {
    lims <- list(c(min(raw$x), max(raw$x)), c(min(raw$y), max(raw$y)))
  }
  x1_lims <- lims[[1]]
  x2_lims <- lims[[2]]
  n2 <- max(10, n2)
  x_seq <- seq(x1_lims[1], x1_lims[2], length = n2)
  y_seq <- seq(x2_lims[1], x2_lims[2], length = n2)
  x_rep <- rep(x_seq, n2)
  y_rep <- rep(y_seq, rep(n2, n2))

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

  list(
    X = X,
    aux = list(
      x = x_seq,
      y = y_seq,
      raw = raw
    )
  )
}


# Internal function for preparing plot of two dimensional smooths
.get_plot_predict_matrix_and_aux_plot_smooth_md <- function(
    mgcv_term,
    fix,
    data = NULL,
    n2 = 40,
    lims = NULL,
    too_far = 0.1,
    ...) {
  ov <- names(fix)
  iv <- mgcv_term$term[!(mgcv_term$term %in% ov)]
  xterm <- iv[1]
  yterm <- iv[2] # TODO use x1, x2 naming convension?
  raw <- data.frame(
    x = as.numeric(data[xterm][[1]]),
    y = as.numeric(data[yterm][[1]])
  )
  n2 <- max(10, n2)

  if (is.null(lims)) {
    lims <- list(c(min(raw$x), max(raw$x)), c(min(raw$y), max(raw$y)))
  }
  x1_lims <- lims[[1]]
  x2_lims <- lims[[2]]
  x1_seq <- seq(x1_lims[1], x1_lims[2], length = n2)
  x2_seq <- seq(x2_lims[1], x2_lims[2], length = n2)

  xx <- rep(x1_seq, n2)
  yy <- rep(x2_seq, each = n2)

  if (mgcv_term$by != "NA") {
    # deal with any by variables
    by <- rep(1, n2^2)
    dat <- data.frame(x = xx, y = yy, by = by) # TODO y is misleading
    colnames(dat) <- c(xterm, yterm, mgcv_term$by)
  } else {
    dat <- data.frame(x = xx, y = yy)
    colnames(dat) <- c(xterm, yterm)
  } ## prediction data.frame complete

  for (ii in ov) {
    dat[[ii]] <- rep(fix[ii], n2^2)
  }

  X <- PredictMat(mgcv_term, dat) ## prediction matrix for this term

  list(
    X = X,
    aux = list(
      x = x1_seq,
      y = x2_seq,
      raw = raw
    )
  )
}
