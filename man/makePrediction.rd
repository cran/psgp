\name{makePrediction}
\alias{makePrediction}
\alias{spatialPredict.psgp}
\title{Spatial projected sequential GP prediction}
\description{\code{makePrediction} is a method for prediction of a variable within the INTAMAP package.
}
\usage{
makePrediction(object, vario)
\method{spatialPredict}{psgp}(object, ...)
}
\arguments{
  \item{object}{ a list object of Intamap type. Most arguments necessary for interpolation
  are passed through this object. 
  See the introduction to the 
  \code{\link[intamap:intamap-package]{intamap-package}} for further 
  description of the necessary content of this variable.  Additional meta data about the measurement process is included in this object.}
  \item{vario}{ an integer array object to describe the variogram to be used for prediction. 
   The first element is a code for the variogram model: 1 - Gaussian, 2 - Exponential.  
   The second element is the range of the process.  
   The third element is the sill.  The fouth element is the noise.  
   The fifth element encodes a constant trend parameter.
  } 
  \item{...}{other parameters to be passed through the generic 
     \code{\link[intamap:spatialPredict]{spatialPredict}} function of the 
      \code{\link[intamap:intamap-package]{intamap-package}}. Possible parameters
       in this function is \code{nsim}, the number of simulations requested.}
} 

\details{
The function \code{makePrediction} is a function for making spatial interpolation with 
Projected Spatial Gaussian Process (PSGP) methods (Csato and Opper, 2002; Ingram et al., 2008).
These methods are able to also take the measurement characteristics into account,
in this function implemented as the element \code{obsChar} in \code{object}.
The parameters can be estimated in \code{\link{learnParameters}}.

Instead of calling this function directly, a user is advised to call the generic S3-class 
wrapper function \code{\link[intamap:spatialPredict]{spatialPredict}} of the
\code{\link[intamap:intamap-package]{intamap-package}} with an \code{object} of class \code{psgp}.


Most of the method is implemented in C++, relying on the external library IT++
(\url{http://itpp.sourceforge.net}), which is a C++ library composed of
classes and functions for linear algebra (matrices and vectors). 

When no measurement metadata is available, the C++ code will default to a Gaussian noise model.

If simulations are requested, these are created with the the function 
\code{\link[gstat:krige]{krige}}-function in the \code{gstat}-package.
}

\references{ 
\url{http://www.intamap.org/}

L. Csato and M. Opper. Sparse online Gaussian processes. Neural Computation, 14(3):
641-669, 2002.

B. Ingram, D. Cornford, and D. Evans. Fast algorithms for automatic mapping with space-
limited covariance functions. Stochastic Environmental Research and Risk Assessment, 22
(5):661-670, 2008.
}
\author{Ben Ingram}
\seealso{
\code{\link[intamap]{spatialPredict}}, \code{\link{learnParameters}}
}
\examples{
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
obj = estimateParameters(obj)
obj = spatialPredict(obj) # directly calls makePrediction using variogram stored in obj
}
\keyword{spatial}
