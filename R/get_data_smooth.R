#' @description Get data for plotting one dimensional smooth effects.
#' @export
get_data.mgcv.smooth.1D <- function(
    term,
    fitted_terms,
    gam,
    n = 100,
    lims = NULL,
    maxpo = 1e4,
    unconditional = FALSE,
    se_with_mean = FALSE,
    nsim = 0,
    ...) {

  P <- prepareP(
    term = term,
    fitted_terms = fitted_terms,
    gam = gam,
    unconditional = unconditional,
    n = n,
    lims = lims,
    se_with_mean = se_with_mean,
    nsim = nsim
  )
  # Construct data.frame of partial residuals
  mgcv_term <- gam$smooth[[term$term_idx]]
  res <- data.frame(x = as.vector(gam$model[mgcv_term$term][[1]]))
  if (!is.null(P$partial_resids) && length(P$partial_resids)) {
    res$y <- P$partial_resids
  }

  fit <- data.frame(
    x = P$aux$x,
    y = P$fit,
    se = P$se
  )
  sim <- NULL
  if (!is.null(P$sim_f)) {
    nsim <- ncol(P$sim_f)
    sim <- data.frame(
      "x" = rep(P$aux$x, nsim),
      "id" = as.factor(rep(1:nsim, each = nrow(P$sim_f)))
    )
  }

  list(
    fit=fit,
    res=res,
    sim=sim,
  )
}

#' @description Get data for plotting two dimensional smooth effects.
#' @export
get_data.mgcv.smooth.2D <- function(
    term,
    gam,
    fitted_terms,
    n = 40,
    lims = NULL,
    se_with_mean = FALSE,
    unconditional = FALSE,
    ...) {

  P <- prepareP(
    term = term,
    fitted_terms = fitted_terms,
    gam = gam,
    unconditional = unconditional,
    n = n,
    lims = lims,
    se_with_mean = se_with_mean
  )
  .get_data_shared_2d(P = P)
}


# Used by e.g MD and SOS
#' @noRd
.get_data_shared_2d <- function(P, flip = FALSE) {
  # 1) Build dataset on fitted effect
  P$fit[P$aux$exclude] <- NA
  fit <- data.frame(
    y = drop(P$fit),
    x1 = rep(P$aux$x1, length(P$fit) / length(P$aux$x1)),
    x2 = rep(P$aux$x2, each = length(P$fit) / length(P$aux$x)),
    se = P$se
  )
  res <- data.frame(
    x1 = P$x_raw$x1,
    x2 = P$x_raw$x2,
    y = P$partial_resids
  )
  if (flip) {
    fit2 <- fit
    res2 <- res
    fit2$x1 <- fit$x2
    fit2$x2 <- fit$x1
    res2$x1 <- res$x2
    res2$x2 <- res$x1
    res <- res2
    fit <- fit2
  }

  list(
    fit = fit,
    res = res
  )
}


#' @description Get data for plotting a 2D slice of a higher-dimensional smooth effects.
#' @export
#'
get_data.mgcv.smooth.MD <- function(
    term,
    gam,
    fitted_terms,
    fix,
    n = 40,
    lims = NULL,
    se_with_mean = FALSE,
    unconditional = FALSE,
    ...) {

  P <- prepareP(
    term = term,
    fitted_terms = fitted_terms,
    gam = gam,
    unconditional = unconditional,
    n = n,
    lims = lims,
    se_with_mean = se_with_mean,
    fix = fix
  )

  .get_data_shared_2d(P = P)
}



#' @description Default plot preparation method for smooth objects `x'
#' inheriting from "mgcv.smooth".
#' @export
.get_plot_predict_matrix_and_aux.mgcv.smooth <- function(
    mgcv_term,
    data = NULL,
    n = NULL,
    lims = NULL,
    ...) {
  if (is.null(n) && mgcv.term$dim==1) {
    n <- 100
  } else if (is.null(n)) {
    n <- 40 # >=2D
  }
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
      n = n,
      lims = lims,
      ...
    )
  }

  if (mgcv_term$dim > 2) {
    out <- .get_plot_predict_matrix_and_aux_plot_smooth_md(
      mgcv_term = mgcv_term,
      data = data,
      n = n,
      lims = lims,
      ...
    )
  }
  out
}


# Internal function for preparing plot of one dimensional smooths
.get_plot_predict_matrix_and_aux_plot_smooth_1d <- function(
    mgcv_term,
    data,
    lims,
    n,
    ...) {
  x_raw <- as.vector(data[mgcv_term$term][[1]])
  x_seq <- seq(lims[1], lims[2], length = n)

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
    aux = list(x = x_seq)
  )
  return(out)
}

# Internal function for preparing plot of two dimensional smooths
.get_plot_predict_matrix_and_aux_plot_smooth_2d <- function(
    mgcv_term,
    lims,
    data = NULL,
    n,
    ...) {
  # 2D case is just MD case with no fixed variables
  .get_plot_predict_matrix_and_aux_plot_smooth_md(
    mgcv_term = mgcv_term,
    fix = NULL,
    data = data,
    n = n,
    lims = lims,
    ...
  )
}


# Internal function for preparing plot of multi-dimensional smooths
.get_plot_predict_matrix_and_aux_plot_smooth_md <- function(
    mgcv_term,
    lims,
    n,
    fix = NULL,
    data = NULL,
    ...) {

  # Determine which terms to vary (exclude fixed variables)
  if (is.null(fix)) {
    # 2D case: use first two terms directly
    x1_term <- mgcv_term$term[1]
    x2_term <- mgcv_term$term[2]
    ov <- character(0)
  } else {
    # MD case: exclude fixed variables
    ov <- names(fix)
    iv <- mgcv_term$term[!(mgcv_term$term %in% ov)] # TODO should we check for length==2?
    x1_term <- iv[1]
    x2_term <- iv[2]
  }

  x_raw <- data.frame(
    x1 = as.numeric(data[x1_term][[1]]),
    x2 = as.numeric(data[x2_term][[1]])
  )
  n <- max(10, n) # TODO this seems likely to produce confusing behaviour, error if >10 provided?

  x1_lims <- lims[[1]]
  x2_lims <- lims[[2]]
  x1_seq <- seq(x1_lims[1], x1_lims[2], length = n)
  x2_seq <- seq(x2_lims[1], x2_lims[2], length = n)

  x1_rep <- rep(x1_seq, n)
  x2_rep <- rep(x2_seq, each = n)

  if (mgcv_term$by != "NA") {
    # deal with any by variables
    by <- rep(1, n^2)
    dat <- data.frame(x = x1_rep, x2 = x2_rep, by = by) # TODO y is misleading
    colnames(dat) <- c(x1_term, x2_term, mgcv_term$by)
  } else {
    dat <- data.frame(x1 = x1_rep, x2 = x2_rep)
    colnames(dat) <- c(x1_term, x2_term)
  } ## prediction data.frame complete

  for (ii in ov) {
    dat[[ii]] <- rep(fix[ii], n^2)
  }

  X <- PredictMat(mgcv_term, dat) ## prediction matrix for this term

  list(
    X = X,
    aux = list(
      x1 = x1_seq,
      x2 = x2_seq,
      x_raw = x_raw
    )
  )
}
