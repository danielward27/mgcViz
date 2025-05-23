#'
#' Extracting parametric effects from a GAM model
#'
#' @description This function can be used to extract a parametric effect from an object of
#'              class \code{gamViz}.
#'
#' @param o an object of class \code{gamViz}, the output of a \code{getViz()} call.
#' @param select index of the selected parametric effect.
#' @return An object of class "pTermSomething" where "Something" is substituted with
#'         the class of the variable of interest. For instance if this "numeric", the \code{pterm}
#'         will return an object of class "ptermNumeric".
#' @name pterm
#' @examples
#' ####### 1. Gaussian GAM
#' library(mgcViz)
#' set.seed(3)
#' dat <- gamSim(1, n = 1500, dist = "normal", scale = 20)
#' dat$fac <- as.factor(sample(c("A1", "A2", "A3"), nrow(dat), replace = TRUE))
#' dat$logi <- as.logical(sample(c(TRUE, FALSE), nrow(dat), replace = TRUE))
#' bs <- "cr"
#' k <- 12
#' b <- gam(y ~ x0 + x1 + I(x1^2) + s(x2, bs = bs, k = k) + fac + x3:fac + I(x1 * x2) + logi, data = dat)
#' o <- getViz(b)
#'
#' # Plot effect of 'x0'
#' pt <- pterm(o, 1)
#' plot(pt, n = 60) + l_ciPoly() + l_fitLine() + l_ciLine() + l_points()
#'
#' # Plot effect of 'x3'
#' pt <- pterm(o, 1)
#' plot(pt, n = 60) + l_fitLine() + l_ciLine(colour = 2)
#'
#' # Plot effect of 'fac'
#' pt <- pterm(o, 4)
#' plot(pt) + l_ciBar(colour = "blue") + l_fitPoints(colour = "red") +
#'   l_rug(alpha = 0.3)
#'
#' # Plot effect of 'logi'
#' pt <- pterm(o, 6)
#' plot(pt) + l_fitBar(a.aes = list(fill = I("light blue"))) + l_ciBar(colour = "blue")
#'
#' # Plot effect of 'x3:fac': no method available yet available for second order terms
#' pt <- pterm(o, 7)
#' plot(pt)
#'
#' ####### 1. Continued: Quantile GAMs
#' b <- mqgamV(y ~ x0 + x1 + I(x1^2) + s(x2, bs = bs, k = k) + x3:fac +
#'   I(x1 * x2) + logi, data = dat, qu = c(0.3, 0.5, 0.8))
#'
#' plot(pterm(b, 3)) + l_ciBar(colour = 2) + l_fitPoints()
#'
#' plot(pterm(b, 4)) + l_fitBar(colour = "blue", fill = 3) + l_ciBar(colour = 2)
#'
#' # Don't know how to plot this interaction
#' plot(pterm(b, 6))
#'
#' ####### 2. Gaussian GAMLSS model
#' library(MASS)
#' mcycle$fac <- as.factor(sample(c("z", "k", "a", "f"), nrow(mcycle), replace = TRUE))
#' b <- gam(list(accel ~ times + I(times^2) + s(times, k = 10), ~ times + fac + s(times)),
#'   data = mcycle, family = gaulss(), optimizer = "efs"
#' )
#' o <- getViz(b)
#'
#' # Plot effect of 'I(times^2)' on mean: notice that partial residuals
#' # are unavailable for GAMLSS models, hence l_point does not do anything here.
#' pt <- pterm(o, 2)
#' plot(pt) + l_ciPoly() + l_fitLine() + l_ciLine() + l_points()
#'
#' # Plot effect of 'times' in second linear predictor.
#' # Notice that partial residuals are unavailable.
#' pt <- pterm(o, 3)
#' plot(pt) + l_ciPoly() + l_fitLine() + l_ciLine(linetype = 3) + l_rug()
#'
#' # Plot effect of 'fac' in second linear predictor.
#' pt <- pterm(o, 4)
#' plot(pt) + l_ciBar(colour = "blue") + l_fitPoints(colour = "red") +
#'   l_rug()
#'
#' @rdname pterm
#' @export pterm
#'
pterm <- function(o, select) {
  if (inherits(o, "list")) {
    if (
      all(sapply(o, function(.x) {
        inherits(.x, "gamViz")
      })) ==
        FALSE
    ) {
      stop(
        "Object \"o\" should be of class \"mgamViz\" or (a list of) \"gamViz\" objects"
      )
    }
    if (is.null(names(o))) {
      names(o) <- 1:length(o)
    }
    class(o) <- "mgamViz"
  }

  if (inherits(o, "mgamViz")) {
    out <- lapply(o, pterm, select = select)
    class(out) <- paste0("multi.", class(out[[1]]))
    attr(out, "isMQGAM") <- inherits(o, "mqgamViz") # Signal that 'o' is output of mqgamV
    return(out)
  }

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
    "ism" = select,
    "name" = nam,
    "varName" = vNam,
    "class" = cls,
    "order" = ord,
    "gObj" = o
  )

  cl <- paste0("pterm", .simpleCap(.mapVarClass(cls)))

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

# Simple function that capitalizes first letter of each word
.simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "", collapse = " ")
}
