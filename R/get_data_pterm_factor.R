#' @description Get the data for plotting of logical parametric effects
#' @export
get_data.pterm_factor <- function(term, fitted_terms, gam, trans = identity, ...) {
  if (term$order > 1) {
    message("mgcViz does not know how to plot this effect. Returning NULL.")
    return(invisible(NULL))
  }

  # 1) Do prediction
  X <- gam$model

  vr <- as.factor(X[[term$varNam]])
  xx <- as.factor(levels(vr))
  data <- X[1:length(xx), ]
  data[[term$varNam]] <- xx

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
    "x" = xx,
    "y" = unname(.pred$fit),
    "ty" = trans(unname(.pred$fit)),
    "se" = unname(.pred$se)
  )
  data$misc <- list("trans" = trans)

  # 3) Get partial residuals
  data$res <- data.frame("x" = as.factor(gam$model[[term$varNam]]))

  # Check if partial residuals are defined: for instance the are not for gamlss models
  if (is.null(gam$residuals) || is.null(gam$weights)) {
    data$res$y <- NULL
  } else {
    .wr <- sqrt(gam$weights)
    .wr <- gam$residuals * .wr / mean(.wr) # weighted working residuals
    data$res$y <- trans(
      .wr +
        fitted_terms[, which(colnames(fitted_terms) == term$nam)]
    )
  }
  data
}
