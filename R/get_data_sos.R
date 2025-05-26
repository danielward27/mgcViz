#' Plotting smooths on the sphere
#'
#' @description This is the plotting method for smooth effects on the sphere.
#' @name plot.sos.smooth
#' @param x a smooth effect object, extracted using [mgcViz::sm].
#' @param scheme if 0 the smooth effect is plotted on the sphere. If 1 the smooth effect is plotted
#'               on the two hemispheres.
#' @param n sqrt of the number of grid points used to compute the effect plot.
#' @param xlim if supplied then this pair of numbers are used as the x limits for the plot.
#' @param ylim if supplied then this pair of numbers are used as the y limits for the plot.
#' @param maxpo maximum number of residuals points that will be used by layers such as
#'              \code{resRug()} and \code{resPoints()}. If number of datapoints > \code{maxpo},
#'              then a subsample of \code{maxpo} points will be taken.
#' @param too_far if greater than 0 then this is used to determine when a location is too far
#'               from data to be plotted. This is useful since smooths tend to go wild
#'               away from data. The data are scaled into the unit square before deciding
#'               what to exclude, and too_far is a distance within the unit square.
#'               Setting to zero can make plotting faster for large datasets, but care
#'               then needed with interpretation of plots.
#' @param phi one of the plotting angles, relevant only if \code{scheme = 0}.
#' @param theta the other plotting angle, relevant only if \code{scheme = 0}.
#' @param trans monotonic function to apply to the smooth and residuals, before plotting.
#'              Monotonicity is not checked.
#' @param unconditional if \code{TRUE} then the smoothing parameter uncertainty corrected covariance
#'                      matrix is used to compute uncertainty bands, if available.
#'                      Otherwise the bands treat the smoothing parameters as fixed.
#' @param seWithMean if TRUE the component smooths are shown with confidence intervals that
#'                   include the uncertainty about the overall mean. If FALSE then the uncertainty
#'                   relates purely to the centred smooth itself. Marra and Wood (2012) suggests
#'                   that TRUE results in better coverage performance, and this is also suggested
#'                   by simulation.
#' @param ... currently unused.
#' @return An objects of class \code{plotSmooth}.
#' @references Marra, G and S.N. Wood (2012) Coverage Properties of Confidence Intervals for
#'             Generalized Additive Model Components. Scandinavian Journal of Statistics.
#' @examples
#' library(mgcViz)
#' set.seed(0)
#' n <- 400
#'
#' f <- function(la, lo) { ## a test function...
#'   sin(lo) * cos(la - .3)
#' }
#'
#' ## generate with uniform density on sphere...
#' lo <- runif(n) * 2 * pi - pi ## longitude
#' la <- runif(3 * n) * pi - pi / 2
#' ind <- runif(3 * n) <= cos(la)
#' la <- la[ind]
#' la <- la[1:n]
#'
#' ff <- f(la, lo)
#' y <- ff + rnorm(n) * .2 ## test data
#'
#' ## generate data for plotting truth...
#' lam <- seq(-pi / 2, pi / 2, length = 30)
#' lom <- seq(-pi, pi, length = 60)
#' gr <- expand.grid(la = lam, lo = lom)
#' fz <- f(gr$la, gr$lo)
#' zm <- matrix(fz, 30, 60)
#'
#' require(mgcv)
#' dat <- data.frame(la = la * 180 / pi, lo = lo * 180 / pi, y = y)
#'
#' ## fit spline on sphere model...
#' bp <- gam(y ~ s(la, lo, bs = "sos", k = 60), data = dat)
#' bp <- getViz(bp)
#'
#' # Plot on sphere
#' plot(sm(bp, 1), scheme = 0) + l_fitRaster() + l_fitContour() +
#'   l_points(shape = 19) + l_rug() + l_coordContour() + l_bound()
#'
#' # Plotting as in standard 2D plots
#' plot(sm(bp, 1), scheme = 1) + l_fitRaster() + l_fitContour() +
#'   l_points(shape = 19) + l_rug()
#' @rdname plot.sos.smooth
#' @export get_data.sos.smooth
#' @export
#'
get_data.sos.smooth <- function(
    term,
    n = 40,
    xlim = NULL,
    ylim = NULL,
    maxpo = 1e4,
    too_far = 0.1,
    phi = 30,
    theta = 30,
    trans = identity,
    scheme = 0,
    seWithMean = FALSE,
    unconditional = FALSE,
    ...) {
  if (length(scheme) > 1) {
    scheme <- scheme[1]
    warning("scheme should be a single number")
  }

  if ((!is.null(xlim) || !is.null(ylim)) && scheme == 0) {
    stop("xlim and ylim must be left to NULL when scheme==0")
  }

  # 1) Prepare data
  P <- prepareP(
    term = term,
    unconditional = unconditional,
    residuals = TRUE,
    resDen = "none",
    se = TRUE,
    se_mult = 1,
    n = NULL,
    n2 = n,
    ylim = ylim,
    xlim = xlim,
    too_far = too_far,
    seWithMean = seWithMean,
    scheme = scheme,
    phi = phi,
    theta = theta
  )

  # 2) Produce output object
  if (scheme == 0) {
    # plot on sphere

    out <- .plot.sos.smooth(term = P$smooth, P = P, trans = trans, maxpo = maxpo)
    out$type <- "sos0"
  } else {
    # standard 2D plot

    out <- .get_data_shared_2d(
      term = P$smooth,
      P = P,
      trans = trans,
      maxpo = maxpo,
      flip = TRUE
    )
    out$type <- "sos1"
  }
  return(out)
}


###########################
# Internal function
.plot.sos.smooth <- function(term, P, trans, maxpo) {
  .dat <- list()

  ### 1) Build dataset on fitted effect
  # We set to NA the entries that fall outside the globe.
  m <- length(P$xm)
  zz <- lo <- la <- se <- rep(NA, m * m)

  se[P$ind] <- P$se
  zz[P$ind] <- P$fit
  lo[P$ind] <- P$lo
  la[P$ind] <- P$la

  .dat$fit <- data.frame(
    "x" = rep(P$xm, m),
    "y" = rep(P$ym, each = m),
    "z" = zz,
    "tz" = trans(zz),
    "lo" = lo,
    "la" = la,
    "se" = se
  )

  ### 2) Build dataset on residuals
  if (!is.null(P$raw)) {
    .dat$res <- P$raw

    # Sample if too many points (> maxpo)
    nres <- nrow(.dat$res)
    .dat$res$sub <- if (nres > maxpo) {
      sample(c(rep(T, maxpo), rep(F, nres - maxpo)))
    } else {
      rep(T, nres)
    }
  }

  .dat$misc <- list("trans" = trans)
  return("data" = .dat)
}


#' @noRd
#' @export
.prepare.sos.smooth <- function(
    term,
    data,
    se1_mult = 1,
    se2_mult = 1,
    partial_resids = NULL,
    se,
    n,
    n2,
    ylim = NULL,
    xlim = NULL,
    too_far,
    trans,
    phi,
    theta,
    scheme,
    ...) {
  ## plot method function for sos.smooth terms
  if (scheme == 1) {

    return(.preparePlotSmooth2D(
      term = term,
      data = data,
      se_mult = se2_mult,
      n2 = n2,
      ylim = ylim,
      xlim = xlim,
      too_far = too_far,
    ))
  }

  ## convert location of pole in plotting grid to radians
  phi <- phi * pi / 180
  theta <- theta * pi / 180

  ## re-map to sensible values...
  theta <- theta %% (2 * pi)
  if (theta > pi) theta <- theta - 2 * pi

  phi <- phi %% (2 * pi)
  if (phi > pi) phi <- phi - 2 * pi
  if (phi > pi / 2) phi <- pi - phi
  if (phi < -pi / 2) phi <- -(phi + pi)

  if (!term$plot.me) {
    return(NULL)
  } ## shouldn't or can't plot
  ## get basic plot data
  raw <- data[term$term]
  raw <- as.data.frame(.lolaxy(
    lo = raw[[2]] * pi / 180,
    la = raw[[1]] * pi / 180,
    theta,
    phi
  ))

  m <- round(n2 * 1.5)
  ym <- xm <- seq(-1, 1, length = m)
  gr <- expand.grid(x = xm, y = ym)
  r <- z <- gr$x^2 + gr$y^2
  z[z > 1] <- NA
  z <- sqrt(1 - z)

  ## generate la, lo in plotting grid co-ordinates...
  ind <- !is.na(z)
  r <- r[ind]
  la <- asin(gr$y[ind])
  lo <- cos(la)
  lo <- asin(gr$x[ind] / lo)

  um <- .repole(lo, la, theta, phi)

  dat <- data.frame(la = um$la * 180 / pi, lo = um$lo * 180 / pi)
  names(dat) <- term$term
  if (term$by != "NA") dat[[term$by]] <- la * 0 + 1

  X <- PredictMat(term, dat) # prediction matrix for this term

  ## fix lo for smooth contouring
  lo <- dat[[2]]
  ii <- lo <= -177
  lo[ii] <- lo[ii] <- 360 + lo[ii]
  ii <- lo < -165 & lo > -177
  ii <- ii | (abs(dat[[1]]) > 80)
  lo[ii] <- NA

  return(list(
    X = X,
    scale = FALSE,
    se = TRUE,
    raw = raw,
    ind = ind,
    xm = xm,
    ym = ym,
    lo = lo,
    la = dat[[1]],
    se_mult = se1_mult
  ))
} ## end prepare.sos.smooth


.repole <- function(lo, la, lop, lap) {
  ## painfully plodding function to get new lo, la relative to pole at
  ## lap,lop...
  ## x,y,z location of pole...
  yp <- sin(lap)
  xp <- cos(lap) * sin(lop)
  zp <- cos(lap) * cos(lop)

  ## x,y,z location of meridian point for pole - i.e. point lat pi/2
  ## from pole on pole's lon.

  ym <- sin(lap - pi / 2)
  xm <- cos(lap - pi / 2) * sin(lop)
  zm <- cos(lap - pi / 2) * cos(lop)

  ## x,y,z locations of points in la, lo

  y <- sin(la)
  x <- cos(la) * sin(lo)
  z <- cos(la) * cos(lo)

  ## get angle between points and new equatorial plane (i.e. plane orthogonal to pole)
  d <- sqrt((x - xp)^2 + (y - yp)^2 + (z - zp)^2) ## distance from points to to pole
  phi <- pi / 2 - 2 * asin(d / 2)

  ## location of images of la,lo on (new) equatorial plane
  ## sin(phi) gives distance to plane, -(xp, yp, zp) is
  ## direction...
  x <- x - xp * sin(phi)
  y <- y - yp * sin(phi)
  z <- z - zp * sin(phi)

  ## get distances to meridian point
  d <- sqrt((x - xm)^2 + (y - ym)^2 + (z - zm)^2)
  ## angles to meridian plane (i.e. plane containing origin, meridian point and pole)...
  theta <- (1 + cos(phi)^2 - d^2) / (2 * cos(phi))
  theta[theta < -1] <- -1
  theta[theta > 1] <- 1
  theta <- acos(theta)

  ## now decide which side of meridian plane...

  ## get points at extremes of hemispheres on either side
  ## of meridian plane....
  y1 <- 0
  x1 <- sin(lop + pi / 2)
  z1 <- cos(lop + pi / 2)

  y0 <- 0
  x0 <- sin(lop - pi / 2)
  z0 <- cos(lop - pi / 2)

  d1 <- sqrt((x - x1)^2 + (y - y1)^2 + (z - z1)^2)
  d0 <- sqrt((x - x0)^2 + (y - y0)^2 + (z - z0)^2)

  ii <- d0 < d1 ## index -ve lon hemisphere
  theta[ii] <- -theta[ii]

  list(lo = theta, la = phi)
} ## end of repole

.lolaxy <- function(lo, la, theta, phi) {
  ## takes locations lo,la, relative to a pole at lo=theta, la=phi.
  ## theta, phi are expressed relative to plotting co-ordinate system
  ## with pole at top. Convert to x,y in plotting co-ordinates.
  ## all in radians!
  er <- .repole(-lo, la, -pi, phi)
  er$lo <- er$lo - theta
  y <- sin(er$la)
  x <- cos(er$la) * sin(er$lo)
  z <- cos(er$la) * cos(er$lo)
  ind <- z < 0
  list(x = x[ind], y = y[ind])
} ## end of lolaxy
