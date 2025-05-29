#' @description Get data for plotting smooth effects on the sphere.

#' @export
#'
get_data.sos.smooth <- function(
    term,
    fitted_terms,
    gam = gam,
    n = 40,
    lims = NULL,
    too_far = 0.1,
    phi = 30,
    theta = 30,
    trans = identity,
    scheme = 0,
    se_with_mean = FALSE,
    unconditional = FALSE,
    ...) {
  if (length(scheme) > 1) {
    scheme <- scheme[1]
    warning("scheme should be a single number")
  }

  if ((!is.null(lims)) && scheme == 0) {
    stop("lims must be left to NULL when scheme==0")
  }

  # 1) Prepare data
  P <- prepareP(
    term = term,
    fitted_terms = fitted_terms,
    gam = gam,
    unconditional = unconditional,
    residuals = TRUE,
    n = NULL,
    n2 = n,
    lims = lims,
    too_far = too_far,
    se_with_mean = se_with_mean,
    scheme = scheme,
    phi = phi,
    theta = theta
  )

  # 2) Produce output object
  if (scheme == 0) {
    # plot on sphere

    out <- .plot.sos.smooth(term = P$aux$smooth, P = P, trans = trans)
    out$type <- "sos0"
  } else {
    # standard 2D plot

    out <- .get_data_shared_2d(
      P = P,
      trans = trans,
      flip = TRUE
    )
    out$type <- "sos1"
  }
  return(out)
}


###########################
# Internal function
.plot.sos.smooth <- function(term, P, trans) {
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
  if (!is.null(P$aux$raw)) {
    .dat$res <- P$aux$raw
  }

  .dat$misc <- list("trans" = trans)
  return("data" = .dat)
}


#' @noRd
#' @export
.get_plot_predict_matrix_and_aux.sos.smooth <- function(
    mgcv_term,
    data,
    n,
    n2,
    lims = NULL,
    too_far,
    trans,
    phi,
    theta,
    scheme,
    ...) {
  ## plot method function for sos.smooth terms
  if (scheme == 1) {
    return(.get_plot_predict_matrix_and_aux_plot_smooth_2d(
      mgcv_term = mgcv_term,
      data = data,
      n2 = n2,
      lims = lims,
      too_far = too_far
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

  ## get basic plot data
  raw <- data[mgcv_term$term]
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
  names(dat) <- mgcv_term$term
  if (mgcv_term$by != "NA") dat[[mgcv_term$by]] <- la * 0 + 1

  X <- PredictMat(mgcv_term, dat) # prediction matrix for this term

  ## fix lo for smooth contouring
  lo <- dat[[2]]
  ii <- lo <= -177
  lo[ii] <- lo[ii] <- 360 + lo[ii]
  ii <- lo < -165 & lo > -177
  ii <- ii | (abs(dat[[1]]) > 80)
  lo[ii] <- NA

  list(
    X = X,
    aux = list(
      raw = raw,
      ind = ind,
      xm = xm,
      ym = ym,
      lo = lo,
      la = dat[[1]]
    )
  )
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
