#' @description This function converts \code{gam} objects into \code{gamViz} objects.
#' @param o an object of class \code{gam}.
#' @export
gam_to_gam_viz <- function(o, nsim = 0, post = FALSE, newdata, ...) {
  if (!inherits(o, "gam")) {
    stop("\"o\" should be of class \"gam\"")
  }

  if (!inherits(o, "gamViz")) {
    tmp <- o$pterms
    np <- if (is.list(tmp)) {
      length(unlist(lapply(tmp, attr, "order")))
    } else {
      length(attr(tmp, "order"))
    }
    ns <- length(o$smooth)
    terms <- predict(o, type = "terms")

    # predict.bam with discrete = T does not predict parametric terms of order 2. Hence we need to add some empty columns.
    nmis <- (np + ns) - ncol(terms)
    if (nmis) {
      M1 <- if (np - nmis) {
        terms[, 1:(np - nmis)]
      } else {
        c()
      }
      M2 <- matrix(
        0,
        nrow(terms),
        nmis,
        dimnames = list(c(), paste(".fakVar", 1:nmis, sep = ""))
      )
      M3 <- if (ns) {
        terms[, (ncol(terms) - ns + 1):ncol(terms)]
      } else {
        c()
      }
      terms <- cbind(M1, M2, M3)
    }

    o$store <- list("termsFit" = terms, "np" = np)

    class(o) <- c("gamViz", class(o))
  }

  # We try to simulate responses. If an error occurs we report it but do no stop.
  # Most likely error if that o$family does not have any simulation method available.
  # NB: we do not allow to use trans() here, as it might lead to problems with check1D and check2D
  if (nsim > 0) {
    if (post) {
      # Posterior simulations OR ...
      tryCatch(
        o$store$sim <- postSim(
          o,
          nsim = nsim,
          newdata = newdata,
          trans = NULL,
          savePar = FALSE,
          ...
        ),
        error = function(e) {
          message(paste("postSim() failed:", e$message))
        }
      )
    } else {
      # ... parameters fixed at MAP
      tryCatch(
        o$store$sim <- simulate(
          o,
          nsim = nsim,
          newdata = newdata,
          trans = NULL,
          ...
        ),
        error = function(e) {
          message(paste("simulate.gam() failed:", e$message))
        }
      )
    }
  }

  # We store new dataset which will be used for checking in check0D, check1D, check2D ...
  # If object 'o' already contained newdata, we either over-write it or set it to NULL
  if (!missing(newdata)) {
    o$store$newdata <- newdata
  } else {
    if (!is.null(o$store$newdata)) {
      message("getViz: newdata removed from gamViz object")
    }
    o$store$newdata <- NULL
  }

  return(o)
}
