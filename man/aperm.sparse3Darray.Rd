\name{aperm.sparse3Darray}
\alias{aperm.sparse3Darray}
\title{
  Transposition of Sparse Array
}
\description{
  Transpose a sparse three-dimensional array by permuting its dimensions.
}
\usage{
   \method{aperm}{sparse3Darray}(a, perm = NULL, resize = TRUE, \dots)
}
\arguments{
  \item{a}{
    A sparse three-dimensional array (object of class \code{"sparse3Darray"}).
  }
  \item{perm}{
    The subscript permutation vector, a permutation of the integers \code{1:3}.
  }
  \item{resize}{
    Logical value specifying whether the dimensions and dimnames of the array
    should also be adjusted, by permuting them according to the permutation.
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  The function \code{\link[base]{aperm}} is generic. This is the method
  for the class \code{"sparse3Darray"} of sparse three-dimensional
  arrays.
}
\value{
  Another sparse three-dimensional array
  (object of class \code{"sparse3Darray"}).
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{sparse3Darray}},
  \code{\link{tensorSparse}}.
}
\examples{
  M <- sparse3Darray(i=1:4, j=sample(1:4, replace=TRUE),
                     k=c(1,2,1,2), x=1:4, dims=c(5,7,2))
  dim(M)
  P <- aperm(M, c(3,1,2))
  dim(P)
}
\keyword{array}
\keyword{manip}
\concept{sparse}
