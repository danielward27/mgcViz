#' @description Get data required for plotting numerical parametric effects.
#' @export
get_data.pterm_numeric <- function(
    term,
    fitted_terms,
    n = 100,
    lims = NULL,
    trans = identity,
    ...) {
  if (term$order > 1) {
    message("mgcViz does not know how to plot this effect. Returning NULL.")
    return(invisible(NULL))
  }

  gam <- term$gam

  # 1) Do prediction
  X <- gam$model

  if (n > nrow(X)) {
    # Model matrix too short, we make it longer
    X <- X[rep(1:nrow(X), ceiling(n / nrow(X))), ]
  }

  if (is.null(lims)) {
    lims <- range(X[[term$varNam]])
  }

  x_seq <- seq(lims[1], lims[2], length = n)
  data <- X[1:n, ]
  data[[term$varNam]] <- x_seq

  # Suppressing spurious warnings from predict.gam
  .pred <- withCallingHandlers(
    predict.gam(
      gam,
      type = "terms",
      se.fit = TRUE,
      terms = term$nam,
      newdata = data
    ),
    warning = function(w) {
      if (
        is.list(gam$formula) &&
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
  data$res <- data.frame("x" = as.vector(gam$model[[term$varNam]]))

  # Check if partial residuals are defined: for instance the are not for gamlss models
  if (is.null(gam$residuals) || is.null(gam$weights)) {
    data$res$y <- NULL
  } else {
    .wr <- sqrt(gam$weights)
    .wr <- gam$residuals * .wr / mean(.wr) # weighted working residuals
    data$res$y <- trans(
      .wr + fitted_terms[, which(colnames(fitted_terms) == term$nam)]
    )
  }

  data$misc <- list("trans" = trans)
  return(data)
}
