#' @description This is the plotting method for random effects (simple random intercepts).
#' @export
get_data.random.effect <- function(term, trans = identity, ...) {
  P <- prepareP(
    term = term,
    unconditional = FALSE,
    residuals = TRUE,
    resDen = "none",
    se = TRUE,
    se_mult = 1,
    n = 100,
    n2 = NULL,
    ylim = NULL,
    xlim = NULL,
    too_far = NULL,
    seWithMean = FALSE
  )
  dat <- list()

  .n <- length(P$fit)
  dat$fit <- data.frame(
    x = qnorm(ppoints(.n)),
    y = sort(P$fit),
    ty = trans(sort(P$fit))
  )

  dat$misc <- list("trans" = trans)
  return(dat)
}


#' @noRd
#' @export
.prepare.random.effect <- function(
    term,
    data = NULL,
    n = 100,
    ylim = NULL,
    xlim = NULL,
    ...) {
  raw <- data[term$term][[1]]
  p <- term$last.para - term$first.para + 1
  X <- diag(p) # prediction matrix for this term

  return(list(
    X = X,
    scale = FALSE,
    se = FALSE,
    raw = raw
  ))
}
