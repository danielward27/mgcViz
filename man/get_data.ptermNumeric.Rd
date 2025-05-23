% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_data_ptermNumeric.R
\name{plot.ptermNumeric}
\alias{plot.ptermNumeric}
\alias{get_data.ptermNumeric}
\title{Plotting numeric parametric effects}
\usage{
\method{get_data}{ptermNumeric}(term, n = 100, xlim = NULL, maxpo = 10000, trans = identity, ...)
}
\arguments{
\item{n}{number of grid points used to compute main effect and c.i. lines.}

\item{xlim}{if supplied then this pair of numbers are used as the x limits for the plot.}

\item{maxpo}{maximum number of residuals points that will be used by layers such as
\code{resRug()} and \code{resPoints()}. If number of datapoints > \code{maxpo},
then a subsample of \code{maxpo} points will be taken.}

\item{trans}{monotonic function to apply to the fit, confidence intervals and residuals,
before plotting. Monotonicity is not checked.}

\item{...}{currently unused.}

\item{x}{a numerical parametric effect object, extracted using \link{pterm}.}
}
\value{
An object of class \code{plotSmooth}.
}
\description{
This is the plotting method for parametric numerical effects.
}
\examples{
# Simulate data and fit GAM
set.seed(3)
dat <- gamSim(1, n = 2000, dist = "normal", scale = 20)
bs <- "cr"
k <- 12
b <- gam(y ~ x0 + x1 + I(x1^2) + s(x2, bs = bs, k = k) +
  I(x1 * x2) + s(x3, bs = bs), data = dat)
o <- getViz(b, nsim = 0)

# Extract first terms and plot it
pt <- pterm(o, 1)
plot(pt, n = 60) + l_ciPoly() + l_fitLine() + l_ciLine()

# Extract I(x1^2) terms and plot it with partial residuals
pt <- pterm(o, 3)
plot(pt, n = 60) + l_ciPoly() + l_fitLine() + l_ciLine() + l_points()

}
