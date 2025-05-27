#'
#' Extracting parametric effects from a GAM model
#'
#' @description This function can be used to extract a parametric effect from an object of
#'              class \code{gamViz}.
#' @export
pterm <- function(o, select) {
  if (!inherits(o, "gamViz")) {
    stop("Argument 'o' should be of class 'gamViz'. See ?getViz")
  }

  terms <- o$pterms
  if (!is.list(terms)) {
    terms <- list(terms)
  }

  order <- lapply(terms, attr, "order")
  np <- sapply(order, length)
  tot <- sum(np)

  if (length(select) > 1) {
    stop("select should be a scalar")
  }
  if (select > tot) {
    stop(paste(
      "select should be smaller than",
      tot,
      "the number of parametric terms in gamObject"
    ))
  }

  vNam <- unlist(sapply(terms, function(.inp) attr(.inp, "term.labels")))[
    select
  ]
  nam <- if (length(terms) > 1) {
    attr(terms, "term.labels")[select]
  } else {
    attr(terms[[1]], "term.labels")[select]
  }

  ord <- unlist(order)[[select]]
  if (ord > 1) {
    # Dealing with interactions OR ...
    cls <- "interaction"
  } else {
    # ... simple effect
    cls <- unlist(sapply(
      terms,
      function(.inp) unname(attr(.inp, "dataClasses"))[-1]
    ))[select]
  }

  if (cls == "ordered") {
    cls <- "factor"
  } # We treat ordered factors as simple factors

  # Here we deal with with smooth effect built "by hand" (e.g. using bs(x0, degree=1))
  if (grepl("nmatrix", cls, fixed = TRUE)) {
    carrier.name <- function(term) {
      if (length(term) > 1L) carrier.name(term[[2L]]) else as.character(term)
    }

    cls <- "matrixNumeric"
    vNam <- carrier.name(str2expression(nam)[[1]]) # Extract name of covariate x
  }

  out <- list(
    term_idx = select,
    name = nam,
    varName = vNam,
    class = cls,
    order = ord,
    gam_viz = o
  )

  cl <- paste0("pterm_", .mapVarClass(cls))

  class(out) <- cl

  return(out)
}


# Extract class of variable
.mapVarClass <- function(.cl) {
  if ("integer" %in% .cl) {
    return("numeric")
  }
  if ("numeric" %in% .cl) {
    return("numeric")
  }
  if ("logical" %in% .cl) {
    return("logical")
  }
  if ("factor" %in% .cl || "character" %in% .cl) {
    return("factor")
  }
  return(.cl) # Not covered by mgcViz
}
