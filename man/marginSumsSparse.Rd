\name{marginSumsSparse}
\alias{marginSumsSparse}
\title{
  Margin Sums of a Sparse Matrix or Sparse Array
}
\description{
  For a sparse matrix or sparse array, compute the sum of array entries
  for a specified margin or margins.
}
\usage{
marginSumsSparse(X, MARGIN)
}
\arguments{
  \item{X}{
    A matrix, an array,
    a sparse matrix (of class \code{"sparseMatrix"} from the
    \pkg{Matrix} package) or a sparse three-dimensional array
    (of class \code{"sparse3Darray"} from the
    \pkg{spatstat.sparse} package).
  }
  \item{MARGIN}{
    Integer or integer vector specifying the margin or margins.
  }
}
\details{
  This function computes the equivalent of
  \code{\link[base]{apply}(X, MARGIN, sum)}
  for sparse matrices and arrays \code{X}. The argument \code{X} may be
  \itemize{
    \item a matrix
    \item an array of any number of dimensions
    \item a sparse matrix (object inheriting class \code{"sparseMatrix"}
    in the \pkg{Matrix} package)
    \item a sparse three-dimensional array (of class
    \code{"sparse3Darray"} 
    from the \pkg{spatstat.sparse} package).
  }
  In the first two cases, the computation is performed by
  calling   \code{\link[base]{apply}(X, MARGIN, sum)}
  and the result is a vector, matrix or array.
  In the last two cases, the result is a single value, a sparse vector,
  a sparse matrix, or a sparse three-dimensional array.
}
\value{
  A single value, vector, matrix, array,
  sparse vector (class \code{"sparseVector"} in the \pkg{Matrix}
  package), sparse matrix (class \code{"sparseMatrix"} in the
  \pkg{Matrix} package), or sparse three-dimensional array
  (class \code{"sparse3Darray"} 
  from the \pkg{spatstat.sparse} package).
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link[base]{apply}}
}
\examples{
  M <- sparse3Darray(i=1:3, j=c(3,1,2), k=4:2,
                     x=round(runif(3), 2), dims=rep(4, 3))
  marginSumsSparse(M, 1:2)
  marginSumsSparse(M, 1)
  marginSumsSparse(M, integer(0)) # equivalent to sum(M)
}
\keyword{array}
\keyword{algebra}
\keyword{manip}
\concept{sparse}
