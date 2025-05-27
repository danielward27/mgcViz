#' @description Get data for plotting one dimensional smooth effects.
#' @export
get_data.mgcv.smooth.1D <- function(
    term,
    n = 100,
    xlim = NULL,
    maxpo = 1e4,
    trans = identity,
    unconditional = FALSE,
    seWithMean = FALSE,
    nsim = 0,
    ...) {
  P <- prepareP(
    term = term,
    unconditional = unconditional,
    residuals = TRUE,
    resDen = "none",
    se = TRUE,
    se_mult = 1,
    n = n,
    n2 = NULL,
    ylim = NULL,
    xlim = xlim,
    too_far = NULL,
    seWithMean = seWithMean,
    nsim = nsim
  )
  .dat <- list()
  if (!is.null(P$raw)) {
    # Construct data.frame of partial residuals
    res <- data.frame("x" = as.vector(P$raw))
    if (!is.null(P$p.resid) & length(P$p.resid)) {
      res$y <- trans(P$p.resid)
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

  if (!is.null(P$simF)) {
    nsim <- ncol(P$simF)
    .dat$sim <- data.frame(
      "x" = rep(P$x, nsim),
      "ty" = trans(as.vector(P$simF)),
      "id" = as.factor(rep(1:nsim, each = nrow(P$simF)))
    )
  }

  .dat$misc <- list("trans" = trans)
  return(.dat)
}

# Internal function for preparing plot of one dimensional smooths
.prepare_plot_smooth_1d <- function(
    term,
    data,
    se_mult = 1,
    n = 100,
    xlim = NULL,
    ...) {
  out <- NULL
  if (term$plot.me) {
    raw <- as.vector(data[term$term][[1]])
    if (is.null(xlim)) {
      # Generate x sequence for prediction
      x_seq <- seq(min(raw), max(raw), length = n)
    } else {
      x_seq <- seq(xlim[1], xlim[2], length = n)
    }
    if (term$by != "NA") {
      # Deal with any by variables
      by <- rep(1, n)
      dat <- data.frame(x = x_seq, by = by)
      names(dat) <- c(term$term, term$by)
    } else {
      dat <- data.frame(x = x_seq)
      names(dat) <- term$term
    } # Finished preparing prediction data.frame
    X <- PredictMat(term, dat) # prediction matrix for this term

    if (is.null(xlim)) {
      xlim <- range(x_seq)
    }
    out <- list(
      X = X,
      x = x_seq,
      scale = TRUE,
      se = TRUE,
      raw = raw,
      se_mult = se_mult,
      xlim = xlim
    )
  }
  return(out)
}

#' @description Default plot preparation method for smooth objects `x'
#' inheriting from "mgcv.smooth".
#' @export
.prepare.mgcv.smooth <- function(
    term,
    data = NULL,
    se1_mult = 1,
    se2_mult = 2,
    n = 100,
    n2 = 40,
    ylim = NULL,
    xlim = NULL,
    too_far = 0.1,
    ...) {
  if (term$dim == 1) {
    out <- .prepare_plot_smooth_1d(
      term = term,
      data = data,
      se_mult = se1_mult,
      n = n,
      xlim = xlim,
      ...
    )
  }

  if (term$dim == 2) {
    out <- .prepare_plot_smooth_2d(
      term = term,
      data = data,
      se_mult = se2_mult,
      n2 = n2,
      ylim = ylim,
      xlim = xlim,
      too_far = too_far,
      ...
    )
  }

  if (term$dim > 2) {
    out <- .prepare_plot_smooth_md(
      term = term,
      data = data,
      se_mult = se2_mult,
      n2 = n2,
      ylim = ylim,
      xlim = xlim,
      too_far = too_far,
      ...
    )
  }

  return(out)
}
