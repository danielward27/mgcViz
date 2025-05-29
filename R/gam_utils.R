#' @export
number_parametric <- function(gam) { # TODO maybe count both term types
  if (is.list(gam$pterms)) {
    n_parametric <- length(unlist(lapply(gam$pterms, attr, "order")))
  } else {
    n_parametric <- length(attr(gam$pterms, "order"))
  }
  n_parametric
}

# TODO Is this added support of newdata fine?
gam_to_fitted_terms <- function(gam, newdata = NULL) {
  n_parametric <- number_parametric(gam)
  n_smooth <- length(gam$smooth)
  if (is.null(newdata)) {
    terms <- predict(gam, type = "terms")
  } else {
    terms <- predict(gam, type = "terms", newdata = newdata)
  }
  # TODO does the bug needing the below workaround still exist?
  # I could not recreate it but perhaps not using correct case.
  # predict.bam with discrete = T does not predict parametric terms
  #  of order 2. Hence we need to add some empty columns.
  n_missing <- (n_parametric + n_smooth) - ncol(terms)
  if (n_missing) {
    M1 <- if (n_parametric - n_missing) {
      terms[, 1:(n_parametric - n_missing)]
    } else {
      c()
    }
    M2 <- matrix(
      0,
      nrow(terms),
      n_missing,
      dimnames = list(c(), paste(".fakVar", 1:n_missing, sep = ""))
    )
    M3 <- if (n_smooth) {
      terms[, (ncol(terms) - n_smooth + 1):ncol(terms)]
    } else {
      c()
    }
    terms <- cbind(M1, M2, M3)
  }
  terms
}
