\name{pitfiller}
\alias{pitfiller}
\alias{pitfiller,SpatRaster-method}
\title{Pit Filler}
\description{
Fills pits (depressions with no outlet) in a digital elevation model (DEM) using the PEM4PIT terrain erosion model.
}
\usage{
\S4method{pitfiller}{SpatRaster}(x, pit = NULL, flowdir = NULL, niter = 10, lambda = 0,
                 deviation_type = "lad", U = 1, D = 300, beta = 0.9,
                 theta_exp = 0.5, filename = "", ...)
}
\arguments{
  \item{x}{SpatRaster with digital elevation model.}
  \item{pit}{SpatRaster with pits. See \code{\link{pitfinder}}.}
  \item{flowdir}{SpatRaster with flow direction or \code{NULL}. If \code{NULL}, it is calculated internally. See \code{\link{flowdirD8lad}}.}
  \item{niter}{Number of iterations. Default is 10.}
  \item{lambda}{Deviation parameter. Default is 0.}
  \item{deviation_type}{Type of deviation. Default is \code{"lad"}. See \code{\link{flowdirD8lad}}.}
  \item{U}{Uplift rate. See PEM4PIT model parameters.}
  \item{D}{Diffusion coefficient. See PEM4PIT model parameters.}
  \item{beta}{Erosion coefficient. See PEM4PIT model parameters.}
  \item{theta_exp}{Exponent for contributing area. See PEM4PIT model parameters.}
  \item{filename}{Character. Output filename.}
  \item{...}{Additional arguments for writing files, as in \code{\link{writeRaster}}.}
}
\value{
A \code{\link{SpatRaster}} object with pits filled.
}
\details{
The PEM4PIT model (Physical Erosion Model for PIT and flat areas correction) corrects DEMs by simulating erosion processes. It improves hydrologic modeling by enhancing flow direction and contributing area representation.

The model is based on the steady-state sediment balance equation:

\deqn{0 = U - \beta A^{\theta} \left( \frac{z - z_d}{\Delta l} \right) + \frac{4D}{\Delta x^2} \left( \bar{z} - z \right)}

Rearranging terms to isolate \eqn{z}:

\deqn{\beta A^{\theta} \left( \frac{z - z_d}{\Delta l} \right) = U + \frac{4D}{\Delta x^2} \left( \bar{z} - z \right)}

Distributing terms:

\deqn{\beta A^{\theta} \frac{z}{\Delta l} - \beta A^{\theta} \frac{z_d}{\Delta l} = U + \frac{4D \bar{z}}{\Delta x^2} - \frac{4D z}{\Delta x^2}}

Combining like terms:

\deqn{\beta A^{\theta} \frac{z}{\Delta l} + \frac{4D z}{\Delta x^2} = U + \frac{4D \bar{z}}{\Delta x^2} + \beta A^{\theta} \frac{z_d}{\Delta l}}

Factoring out \eqn{z}:

\deqn{z \left( \beta A^{\theta} \frac{1}{\Delta l} + \frac{4D}{\Delta x^2} \right) = U + \frac{4D \bar{z}}{\Delta x^2} + \beta A^{\theta} \frac{z_d}{\Delta l}}

Solving for \eqn{z}:

\deqn{z = \frac{U \Delta l \Delta x^2 + 4D \Delta l \bar{z} + \beta A^{\theta} \Delta x^2 z_d}{\beta A^{\theta} \Delta x^2 + 4D \Delta l}}
}
\references{
Grimaldi, S., Nardi, F., Di Benedetto, F., Istanbulluoglu, E., & Bras, R. L. (2007). A physically-based method for removing pits in digital elevation models. \emph{Advances in Water Resources}, 30(10), 2151--2158. \doi{10.1016/j.advwatres.2006.11.016}

Tucker, G. E., & Bras, R. L. (1998). Hillslope processes, drainage density, and landscape morphology. \emph{Water Resources Research}, 34(10), 2751--2764. \doi{10.1029/98WR01474}

Niemann, J. D., Bras, R. L., & Veneziano, D. (2003). A physically based interpolation method for fluvially eroded topography. \emph{Water Resources Research}, 39, 1017. \doi{10.1029/2001WR001050}

Niemann, J. D., Bras, R. L., Veneziano, D., & Rinaldo, A. (2001). Impacts of surface elevation on the growth and scaling properties of simulated river networks. \emph{Geomorphology}, 40(1--2), 37--55. \doi{10.1016/S0169-555X(01)00036-8}
}
\author{
Emanuele Cordano
}
\seealso{
\code{\link{terrain}}, \code{\link{watershed}}, \code{\link{flowdirD8lad}}, \code{\link{pitfinder}}
}
\examples{
f <- system.file("ex/elev.tif", package = "terra")
elev <- rast(f) |> project(y = "epsg:32632")
flowdir <- flowdirD8lad(elev, lambda = 0)
pits <- pitfinder(flowdir, pits_on_boundary = FALSE)
elev2 <- pitfiller(x = elev, pit = pits)
flowdir2 <- terrain(elev2, "flowdir")
flowdir2 <- flowdirD8lad(elev2, lambda = 0)
pits2 <- pitfinder(flowdir, pits_on_boundary = FALSE)
}
\keyword{spatial}
