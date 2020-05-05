\name{tensorSparse}
\alias{tensorSparse}
\title{
  Tensor Product of Sparse Vectors, Matrices or Arrays
}
\description{
  Compute the tensor product of two vectors, matrices or arrays
  which may be sparse or non-sparse. 
}
\usage{
tensorSparse(A, B, alongA = integer(0), alongB = integer(0))
}
\arguments{
  \item{A,B}{
    Vectors, matrices, three-dimensional arrays,
    or objects of class \code{sparseVector}, \code{sparseMatrix}
    or \code{sparse3Darray}.
  }
  \item{alongA}{
    Integer vector specifying the dimensions of \code{A} to be collapsed.
  }
  \item{alongB}{
    Integer vector specifying the dimensions of \code{B} to be collapsed.
  }
}
\details{
  This function is a generalisation, to sparse arrays, of the function
  \code{\link[tensor]{tensor}} in the \pkg{tensor} package.

  \code{tensorSparse} has the same syntax and interpretation
  as \code{\link[tensor]{tensor}}.
  For example, if \code{A} and \code{B} are matrices, then
  \code{tensor(A,B,2,1)} is the matrix product \code{A \%*\% B}
  while \code{tensor(A,B,2,2)} is \code{A \%*\% t(B)}.

  This function \code{tensorSparse} handles
  sparse vectors (class \code{"sparseVector"} in the \pkg{Matrix} package), 
  sparse matrices (class \code{"sparseMatrix"} in the \pkg{Matrix}
  package) and sparse three-dimensional arrays (class
  \code{"sparse3Darray"} in the \pkg{spatstat.sparse} package)
  in addition to the usual vectors, matrices and arrays.

  The result is a sparse object if at least one of \code{A} and \code{B}
  is sparse. Otherwise, if neither \code{A} nor \code{B} is sparse, then the
  result is computed using \code{\link[tensor]{tensor}}.

  The main limitation is that the result cannot have more than 3
  dimensions (because sparse arrays with more than 3 dimensions are not yet
  supported).
}
\value{
  Either a scalar, a vector, a matrix, an array,
  a sparse vector (class \code{"sparseVector"} in the \pkg{Matrix} package), 
  a sparse matrix (class \code{"sparseMatrix"} in the \pkg{Matrix}
  package) or a sparse three-dimensional array (class
  \code{"sparse3Darray"} in the \pkg{spatstat.sparse} package).
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{sparse3Darray}},
  \code{\link{aperm.sparse3Darray}}
}
\examples{
  M <- sparse3Darray(i=1:4, j=sample(1:4, replace=TRUE),
                     k=c(1,2,1,2), x=1:4, dims=c(5,5,2))
  A <- tensorSparse(M, M, 1:2, 2:1)
}
\keyword{array}
\keyword{algebra}
\concept{sparse}

