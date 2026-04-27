\name{symmatrix}
\alias{symmatrix}
\title{
  Create a Symmetric Matrix from a Triangular Subset
}
\description{
  Given a vector containing values for the upper or lower triangular
  entries of a matrix, reconstruct the original matrix, assuming it
  is symmetric.
}
\usage{
symmatrix(x, from = c("lower", "upper"), diag = TRUE, na.value=NA)
}
\arguments{
  \item{x}{
    A vector or sparse vector of atomic values.
  }
  \item{from}{
    Character string (partially matched) specifying whether the
    entries of \code{x} represent the upper or lower triangular
    entries of the matrix.
  }
  \item{diag}{
    Logical value specifying whether \code{x} includes diagonal entries
    of the matrix.
  }
  \item{na.value}{
    Value to be assigned to diagonal entries of the matrix, if
    \code{diag=FALSE}.
  }
}
\details{
  This function reconstructs a symmetric matrix from its
  lower triangular entries, or from its upper triangular entries.
  \itemize{
    \item
    If \code{from="lower"} (the default), the vector \code{x} is assumed
    to contain the lower-triangular entries of a matrix \code{M} (that is, the
    entries \code{M[i,j]} in row \eqn{i} and column \eqn{j} where
    \eqn{i \ge j}{i >= j})
    listed in column-major order (that is, all of column 1 first, then all
    the relevant entries in column 2, etc.). This vector could have been
    obtained from the matrix by \code{x <- M[lower.tri(M, diag)]}.
    \item
    If \code{from="upper"}, the vector \code{x} is assumed
    to contain the upper-triangular entries (those for which
    \eqn{i \le j}{i <= j}) listed in column-major order.
    This vector could have been
    obtained from the matrix by \code{x <- M[upper.tri(M, diag)]}.
  }
  Assuming \code{M} was a
  symmetric matrix, it will be reconstructed from these values.
  The reconstructed matrix is returned.

  If \code{diag} is FALSE, then the diagonal entries of the matrix
  are assumed to have been omitted, and the
  resulting matrix will contain \code{NA} entries in the diagonal.

  If \code{x} is a sparse vector (inheriting class
  \code{"sparseVector"}) the result is a sparse matrix
  (inheriting classes \code{"sparseMatrix"} and \code{"symmetricMatrix"}).
}
\value{
  A symmetric, square matrix or sparse matrix,
  containing values of the same atomic type as \code{x}.
}
\author{
  \adrian.
}
\seealso{
  \code{\link[base]{lower.tri}}
}
\examples{
  M <- matrix(c(10,  2, 3, 4,
                 2, 20, 5, 6,
                 3,  5, 30, 7,
                 4,  6, 7, 40), 4, 4)
  M

  (x <- M[lower.tri(M, diag=TRUE)])
  symmatrix(x)

  xsp <- as(x, "sparseVector")
  symmatrix(xsp)

  (y <- M[upper.tri(M, diag=FALSE)])
  symmatrix(y, from="upper", diag=FALSE)

}
\keyword{manip}
\keyword{array}
