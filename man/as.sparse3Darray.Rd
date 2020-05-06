\name{as.sparse3Darray}
\alias{as.sparse3Darray}
\title{
  Convert Data to a Sparse Three-Dimensional Array
}
\description{
  Convert other kinds of data to a sparse three-dimensional array.
}
\usage{
as.sparse3Darray(x, \dots)
}

\arguments{
  \item{x}{
    Data in another format (see Details).
  }
  \item{\dots}{
    Ignored.
  }
}
\details{
  This function converts data in various formats
  into a sparse three-dimensional array (object of class
  \code{"sparse3Darray"}).

  The argument \code{x} can be
  \itemize{
    \item a sparse three-dimensional array (class \code{"sparse3Darray"})
    \item an \code{array}
    \item a \code{matrix},
    which will be interpreted as an array
    with dimension \code{c(dim(x), 1)}
    \item a sparse matrix (inheriting class \code{"sparseMatrix"} in the
    \pkg{Matrix} package)
    which will be interpreted as an array
    with dimension \code{c(dim(x), 1)}
    \item a vector of atomic values,
    which will be interpreted as an array of dimension
    \code{c(length(x), 1, 1)}
    \item a sparse vector (inheriting class \code{"sparseVector"} in the
    \pkg{Matrix} package)
    which will be interpreted as an array of dimension \code{c(x@length, 1, 1)}
    \item a list of matrices with the same dimensions,
    which will be interpreted as slices \code{A[,,k]} of an array \code{A}
    \item a list of sparse matrices (each inheriting class
    \code{"sparseMatrix"} in the \pkg{Matrix} package) with the same
    dimensions,
    which will be interpreted as slices \code{A[,,k]} of an array \code{A}.
  }
}
\value{
  Sparse three-dimensional array (object of class \code{"sparse3Darray"}).
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{sparse3Darray}}
}
\examples{
  A <- array(c(1,3,0,0,0,0,0,4,0,2,0,5,
               0,0,1,0,0,0,1,0,0,0,1,0),
             dim=c(3,4,2))
  #' array to sparse array
  B <- as.sparse3Darray(A) # positive extent
  #' list of matrices to sparse array
  B <- as.sparse3Darray(list(A[,,1], A[,,2]))
  #' matrix to sparse array
  B1 <- as.sparse3Darray(A[,,1])
  #' vector to sparse array
  B11 <- as.sparse3Darray(A[,1,1])
}
\keyword{array}
\keyword{manip}
\concept{sparse}
