#' @description Get the data for plotting of logical parametric effects
#' @export
get_data.pterm_factor <- function(term, trans = identity, ...) {
  if (term$order > 1) {
    message("mgcViz does not know how to plot this effect. Returning NULL.")
    return(invisible(NULL))
  }

  gam_viz <- term$gam_viz

  # 1) Do prediction
  X <- gam_viz$model

  vr <- as.factor(X[[term$varNam]])
  xx <- as.factor(levels(vr))
  data <- X[1:length(xx), ]
  data[[term$varNam]] <- xx

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
    "x" = xx,
    "y" = unname(.pred$fit),
    "ty" = trans(unname(.pred$fit)),
    "se" = unname(.pred$se)
  )
  data$misc <- list("trans" = trans)

  # 3) Get partial residuals
  data$res <- data.frame("x" = as.factor(gam_viz$model[[term$varNam]]))

  # Check if partial residuals are defined: for instance the are not for gamlss models
  if (is.null(gam_viz$residuals) || is.null(gam_viz$weights)) {
    data$res$y <- NULL
  } else {
    .wr <- sqrt(gam_viz$weights)
    .wr <- gam_viz$residuals * .wr / mean(.wr) # weighted working residuals
    data$res$y <- trans(
      .wr +
        gam_viz$store$termsFit[, which(colnames(gam_viz$store$termsFit) == term$nam)]
    )
  }
  data
}
