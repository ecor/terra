\name{panel}

\docType{methods}

\alias{panel}
\alias{panel,SpatRaster-method}

\title{Make a panel-map}

\description{
Make a panel map from the layers of a SpatRaster
}

\usage{
\S4method{panel}{SpatRaster,missing}(x, main, loc.main="topleft", nc, nr, maxnl=16, maxcell=500000, box=FALSE, pax=list(), plg=list(), ...)
}

\arguments{
  \item{x}{SpatRaster}
  \item{main}{character. Main plot titles (one for each layer to be plotted). You can use arguments \code{cex.main}, \code{font.main}, \code{col.main} to change the appearance}
  \item{loc.main}{numeric of character to set the location of the main title. Either two coordinates, or a character value such as "topleft")}
  \item{nc}{positive integer. Optional. The number of columns to divide the plotting device in (when plotting multiple layers)}
  \item{nr}{positive integer. Optional. The number of rows to divide the plotting device in (when plotting multiple layers)}
  \item{maxnl}{positive integer. Maximum number of layers to plot (for a multi-layer object)} 
  \item{maxcell}{positive integer. Maximum number of cells to use for the plot}
  \item{box}{logical. Should a box be drawn around the map?}
  \item{plg}{see \code{\link{plot}}}
  \item{pax}{see \code{\link{plot}}}
  \item{...}{arguments passed to \code{plot("SpatRaster", "numeric")} and additional graphical arguments}  
}

\seealso{
 \code{\link{plot}} and see \code{rasterVis::levelplot} and \code{tidyterra::autoplot} for more sophisticated panel plots.
}

\examples{ 
f <- system.file("ex/elev.tif", package="terra") 
r <- rast(f)
r <- c(r, r/2, r*2, r)
panel(r)
}


\keyword{methods}
\keyword{spatial}