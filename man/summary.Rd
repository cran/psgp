\name{summary.psgp}
\alias{summary.psgp}
\alias{plot.psgp}
\title{plot and summary functions for psgp}
\description{plot and summary function for psgp
}
\usage{
\method{summary}{psgp}(object, ...)
\method{plot}{psgp}(x, ...)
}
\arguments{
\item{object}{  a list object. Most arguments necessary for interpolation
are passed through this object. See \code{\link[intamap:intamap-package]{intamap-package}} for further
description of the necessary content of this variable}
\item{x}{  same as object}
\item{...}{parameters to be passed to the default summary function for each element, or
   to \code{\link[intamap:plotIntamap]{plotIntamap}}}
}

\details{
These are new methods for \code{psgp}-class objects. The functions call
\code{\link[intamap:plotIntamap]{plotIntamap}} and \code{\link[intamap:summaryIntamap]{summaryIntamap}},
respectively.
}

\references{
\url{http://www.intamap.org/}
}
\seealso{
\code{\link[intamap]{spatialPredict}}, \code{\link{learnParameters}}
}
\keyword{spatial}
