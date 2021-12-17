\name{bind.sparse3Darray}
\alias{bind.sparse3Darray}
\title{
  Combine Three-Dimensional Sparse Arrays
}
\description{
  Two sparse arrays will be joined to make a larger sparse array.
}
\usage{
bind.sparse3Darray(A, B, along)
}
\arguments{
  \item{A,B}{
    Sparse three-dimensional arrays
    (objects of class \code{"sparse3Darray"})
    or data acceptable to \code{\link{as.sparse3Darray}}.
  }
  \item{along}{
    The dimension along which the two arrays will be joined.
    An integer from 1 to 3. 
  }
}
\details{
  This operation is similar to \code{\link[base]{rbind}}, \code{\link[base]{cbind}}
  and \code{\link[abind]{abind}}. The two 3D arrays \code{A} and \code{B}
  will be joined to make a larger 3D array by concatenating them along
  the dimension specified by \code{along}.

  The arguments \code{A} and \code{B} should be sparse three-dimensional arrays
  (objects of class \code{"sparse3Darray"})
  or data acceptable to \code{\link{as.sparse3Darray}}.
  They must have identical array dimensions except in the dimension
  specified by \code{along}.
}
\value{
  A sparse three-dimensional array (object of class \code{"sparse3Darray"}).
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{as.sparse3Darray}},
  \code{\link{methods.sparse3Darray}}.

  See \code{\link[abind]{abind}} for joining non-sparse arrays.
}
\examples{
  M <- sparse3Darray(i=1:3, j=c(3,1,2), k=4:2,
                     x=runif(3), dims=rep(4, 3))
  dim(M)
  U <- M[ , 1:3, ]
  dim(U)
  V <- bind.sparse3Darray(M, U, along=2)
  dim(V)
}
\keyword{array}
\keyword{sparse}
\keyword{manip}
