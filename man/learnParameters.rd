\name{learnParameters}
\alias{learnParameters}
\alias{estimateParameters.psgp}
\title{Projected Sequential Gaussian Process}
\description{\code{learnParameters} performs maximum likelihood parameter estimation in the PSGP framework.
}
\usage{
learnParameters(object)
\method{estimateParameters}{psgp}(object, ...)
}
\arguments{
  \item{object}{ a list object of Intamap type. Most arguments necessary for interpolation
  are passed through this object. See the introduction to the 
  \code{\link[intamap:intamap-package]{intamap-package}} for further 
  description of the necessary content of this variable}
  \item{...}{other parameters possible to pass to \code{\link[intamap:estimateParameters]{estimateParameters}}
      not relevant for \code{estimateParameters.psgp} }
}
 

\details{
The function \code{learnParameters} is a function for estimating variogram parameters with 
Projected Spatial Gaussian Process (PSGP) methods (Csato and Opper, 2002; Ingram et al., 2008)
through a maximum likelihood estimation. 

Instead of calling this function directly, a user is advised to call the generic S3-class 
wrapper function \code{\link[intamap:estimateParameters]{estimateParameters}} of the
\code{\link[intamap:intamap-package]{intamap-package}} with an \code{object} of class \code{psgp}.

Predictions can be done with 
\code{\link{makePrediction}} or \code{\link[intamap:spatialPredict]{spatialPredict}} with
an \code{object} of type \code{psgp}.
These methods are able to also take the measurement characteristics into account,
in this function implemented as the element \code{obsChar} in \code{object}.

Most of the method is implemented in C++, relying on the external library IT++
(\url{http://itpp.sourceforge.net}), which is a C++ library composed of
classes and functions for linear algebra (matrices and vectors). 

}


\references{ 

\url{http://www.intamap.org/}

}
\author{Ben Ingram}
\seealso{
\code{\link[intamap:estimateParameters]{estimateParameters}},\code{\link{makePrediction}}
}
\examples{
# This example skips some steps that might be necessary for more complicated
# tasks, such as estimateParameters and pre- and postProcessing of the data
data(meuse)
coordinates(meuse) = ~x+y
meuse$value = log(meuse$zinc)
data(meuse.grid)
gridded(meuse.grid) = ~x+y
proj4string(meuse) = CRS("+init=epsg:28992")
proj4string(meuse.grid) = CRS("+init=epsg:28992")

# set up intamap object:
obj = createIntamapObject(
	observations = meuse,
	predictionLocations = meuse.grid,
	targetCRS = "+init=epsg:3035",
	class = "psgp"
)

# do interpolation step:
obj = conformProjections(obj)
obj = estimateParameters(obj)  # obj updated with variogram
}
\keyword{spatial}
