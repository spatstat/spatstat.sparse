\name{as.array.sparse3Darray}
\alias{as.array.sparse3Darray}
\title{
  Convert Sparse Array to Full Array
}
\description{
  Convert a sparse three-dimensional array to a full three-dimensional array.
}
\usage{
\method{as.array}{sparse3Darray}(x, \dots)
}
\arguments{
  \item{x}{
    Sparse three-dimensional array
    (object of class \code{"sparse3Darray"}).
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  This is a method for the generic \code{\link[base]{as.array}}
  for sparse three-dimensional arrays (class \code{"sparse3Darray"}).
  It converts the sparse three-dimensional array \code{x} into
  an \code{\link[base]{array}} representing the same data.
}
\value{
  An array (class \code{"array"}) with the same dimensions as \code{x}
  and the same type of entries as \code{x}.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{sparse3Darray}},
  \code{\link{as.sparse3Darray}}
}
\examples{
  M <- sparse3Darray(i=1:3, j=c(3,1,2), k=4:2,
                     x=runif(3), dims=rep(4, 3))
  V <- as.array(M)
}
\keyword{array}
\keyword{manip}
\concept{sparse}
