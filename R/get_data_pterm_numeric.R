#' @description Get data required for plotting numerical parametric effects.
#' @export
get_data.pterm_numeric <- function(
    term,
    n = 100,
    xlim = NULL,
    maxpo = 1e4,
    trans = identity,
    ...) {
  if (term$order > 1) {
    message("mgcViz does not know how to plot this effect. Returning NULL.")
    return(invisible(NULL))
  }

  gam_viz <- term$gam_viz

  # 1) Do prediction
  X <- gam_viz$model

  if (n > nrow(X)) {
    # Model matrix too short, we make it longer
    X <- X[rep(1:nrow(X), ceiling(n / nrow(X))), ]
  }

  if (is.null(xlim)) {
    xlim <- range(X[[term$varNam]])
  }

  x_seq <- seq(xlim[1], xlim[2], length = n)
  data <- X[1:n, ]
  data[[term$varNam]] <- x_seq

  # Suppressing spurious warnings from predict.gam
  .pred <- withCallingHandlers(
    predict.gam(
      gam_viz,
      type = "terms",
      se.fit = TRUE,
      terms = term$nam,
      newdata = data
    ),
    warning = function(w) {
      if (
        is.list(gam_viz$formula) &&
          any(grepl("is absent, its contrast will be ignored", w))
      ) {
        invokeRestart("muffleWarning")
      }
    }
  )

  # 2) Build dataset on fitted effect
  data <- list()
  data$fit <- data.frame(
    "x" = x_seq,
    "y" = unname(.pred$fit),
    "ty" = trans(unname(.pred$fit)),
    "se" = unname(.pred$se)
  )

  # 3) Get partial residuals
  data$res <- data.frame("x" = as.vector(gam_viz$model[[term$varNam]]))

  # Check if partial residuals are defined: for instance the are not for gamlss models
  if (is.null(gam_viz$residuals) || is.null(gam_viz$weights)) {
    data$res$y <- NULL
  } else {
    .wr <- sqrt(gam_viz$weights)
    .wr <- gam_viz$residuals * .wr / mean(.wr) # weighted working residuals
    data$res$y <- trans(
      .wr + gam_viz$store$termsFit[, which(colnames(gam_viz$store$termsFit) == term$nam)]
    )
  }

  # Exclude residuals falling outside boundaries
  data$res <- data$res[
    data$res$x >= xlim[1] & data$res$x <= xlim[2], ,
    drop = FALSE
  ]

  # Sample if too many points (> maxpo)
  nres <- nrow(data$res)
  data$res$sub <- if (nres > maxpo) {
    sample(c(rep(T, maxpo), rep(F, nres - maxpo)))
  } else {
    rep(T, nres)
  }

  data$misc <- list("trans" = trans)
  return(data)
}
