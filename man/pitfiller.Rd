\docType{methods}
\name{pitfiller}

\alias{pitfiller}
\alias{pitfiller,SpatRaster-method}


\title{Pit Filler}

\description{
It fills the pits (depressions with no outlet ) through a terrain erosion model (PEM4PIT)
}

\usage{
\S4method{pitfiller}{SpatRaster}(x, pit = NULL, flowdir = NULL, niter = 10, lambda = 0,
                 deviation_type = "lad", U = 1, D = 300, beta = 0.9,
                 theta_exp = 0.5, filename = "", ...)

}

\arguments{
  \item{x}{SpatRaster wih digital elevation model.} 
  \item{pit}{SpatRaster wih pits. See \code{\link{pitfinder}}}
  \item{flowdir}{SpatRaster wih flow direction or \code{NULL}, intarnally calculated. See \code{\link{flowdirD8lad}}}
  \item{niter}{number of iterations}
  \item{lambda,deviation_type}{default is \code{deviation_type="lad"} and \code{lambda=0}. see \code{\link{flowdirD8lad}} .}  
  \item{U,D,beta,theta_exp}{See PEM4PIT model parameters. See details. } 
  \item{filename}{character. Output filename}
  \item{...}{additional arguments for writing files as in \code{\link{writeRaster}}}
}

\value{ 

A \code{\link{SpatRaster-class}} (raster) delevation model with pits filled. 
}

\details{
  The PEM4PIT model, which stands for Physical Erosion Model for PIT and flat areas correction, is designed to address issues in Digital Elevation Models (DEMs) related to pits and flat areas. These issues can significantly affect hydrologic modeling and analysis. The model uses a simplified physical representation of erosion processes at the basin scale to correct these inaccuracies.

  - **Automatic parameter estimation**: It automatically calibrates parameters to improve the accuracy of DEMs.
  - **Hydrologic applications**: It enhances the representation of flow directions and contributing areas, which are crucial for hydrologic simulations.
  - **Erosion processes**: The model accounts for both fluvial and diffusive erosion processes, providing a more realistic depiction of terrain.

  Model is based on the steady sediment balance equation:
  
  \deqn{0 = U - \beta A^{\theta} \left( \frac{z - z_d}{\Delta l} \right) + \frac{4D}{\Delta x^2} \left( \overset{\smile}{z} - z \right)}

  where:
  
  - \eqn{U}: uplift
  - \eqn{\beta}: A coefficient that represents the influence of the erosion process on the terrain.
  - \eqn{A}: The area contributing to the flow, which affects the erosion rate.
  - \eqn{\theta}: An exponent that modifies the effect of the contributing area \eqn{A} on the erosion process.
  - \eqn{z}: The elevation at a specific point in the DEM.
  - \eqn{z_d}: The elevation of a downstream point, used to calculate the slope.
  - \eqn{\Delta l}: The distance between the points \eqn{z} and \eqn{z_d}, representing the length over which the slope is calculated.
  - \eqn{D}: A diffusion coefficient that represents the rate of terrain smoothing due to diffusive erosion processes.
  - \eqn{\Delta x}: The grid spacing in the DEM, used to normalize the diffusion term.
  - \eqn{\overset{\smile}{z}}: The average elevation of neighboring points, used to calculate the diffusive erosion effect.

  ### Solving the Equation for \eqn{z}

  Starting with the equation:
  
  \deqn{0 = U - \beta A^{\theta} \left( \frac{z - z_d}{\Delta l} \right) + \frac{4D}{\Delta x^2} \left( \overset{\smile}{z} - z \right)}

  1. Rearrange the terms to isolate \eqn{z}:
  
  \deqn{\beta A^{\theta} \left( \frac{z - z_d}{\Delta l} \right) = U + \frac{4D}{\Delta x^2} \left( \overset{\smile}{z} - z \right)}

  2. Distribute \eqn{\beta A^{\theta}} and \eqn{\frac{4D}{\Delta x^2}}:
  
  \deqn{\beta A^{\theta} \frac{z}{\Delta l} - \beta A^{\theta} \frac{z_d}{\Delta l} = U + \frac{4D \overset{\smile}{z}}{\Delta x^2} - \frac{4D z}{\Delta x^2}}

  3. Combine like terms involving \eqn{z}:
  
  \deqn{\beta A^{\theta} \frac{z}{\Delta l} + \frac{4D z}{\Delta x^2} = U + \frac{4D \overset{\smile}{z}}{\Delta x^2} + \beta A^{\theta} \frac{z_d}{\Delta l}}

  4. Factor out \eqn{z}:
  
  \deqn{z \left( \beta A^{\theta} \frac{1}{\Delta l} + \frac{4D}{\Delta x^2} \right) = U + \frac{4D \overset{\smile}{z}}{\Delta x^2} + \beta A^{\theta} \frac{z_d}{\Delta l}}

  5. Solve for \eqn{z}:
  
  \deqn{z = \frac{U \Delta l \Delta x^2 + 4D \Delta l \overset{\smile}{z} + \beta A^{\theta} \Delta x^2 z_d}{\beta A^{\theta} \Delta x^2 + 4D \Delta l}}

  So, the solution for \eqn{z} is:
  
  \deqn{z = \frac{U \Delta l \Delta x^2 + 4D \Delta l \overset{\smile}{z} + \beta A^{\theta} \Delta x^2 z_d}{\beta A^{\theta} \Delta x^2 + 4D \Delta l}}

  see reference: 
  Grimaldi, S., Nardi, F., Di Benedetto, F., Istanbulluoglu, E., & Bras, R. L. (2007). 
  "A physically-based method for removing pits in digital elevation models."
  Advances in Water Resources, 30(10), 2151-2158. doi:10.1016/j.advwatres.2006.11.016
  \url{https://www.sciencedirect.com/science/article/pii/S0309170807000644}
  
  
  Further references: 
  
  Tucker, G. E., and R. L. Bras (1998), Hillslope processes, drainage density, and landscape morphology, Water Resour. Res., 34(10), 2751–2764, doi:10.1029/98WR01474 \url{https://agupubs.onlinelibrary.wiley.com/doi/abs/10.1029/98wr01474} .
  
  Niemann, J. D., R. L. Bras, and D. Veneziano (2003), A physically based interpolation method for fluvially eroded topography, Water Resour. Res., 39, 1017, doi:10.1029/2001WR001050, 1. \url{https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2001WR001050}
  
  Niemann, J. D., R. L. Bras, D. Veneziano and Rinaldo A. (2001) Jeffrey D. Niemann, Rafael L. Bras, Daniele Veneziano, Andrea Rinaldo, Impacts of surface elevation on the growth and scaling properties of simulated river networks,Geomorphology,Volume 40, Issues 1–2,2001,Pages 37-55,ISSN 0169-555X,
https://doi.org/10.1016/S0169-555X(01)00036-8.
\url{https://www.sciencedirect.com/science/article/pii/S0169555X01000368} Keywords: Drainage networks; Models; Self-similarity; Fractal geometry
  
}


\author{

Emanuele Cordano
}

\seealso{\code{\link{terrain}},\code{\link{watershed}},\code{\link{flowdirD8ltd}},\code{\link{pitfinder}}}
\examples{

f <- system.file("ex/elev.tif", package="terra")
elev <- rast(f) |> project(y="epsg:32632")
flowdir <- terrain(elev,"flowdir")
pits <- pitfinder(flowdir)+elev*0

elev2 <- pitfiller(x=elev,pit=pits)


flowdir2 <- terrain(elev2,"flowdir")
pits2 <- pitfinder(flowdir2)







}




\keyword{spatial}
