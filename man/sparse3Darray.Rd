\name{sparse3Darray}
\alias{sparse3Darray}
\title{
  Create a Sparse Three-Dimensional Array
}
\description{
  Create a sparse representation of a three-dimensional array.
}
\usage{
sparse3Darray(i = integer(0), j = integer(0), k = integer(0),
              x = numeric(0),
              dims = c(max(i), max(j), max(k)), dimnames = NULL,
              strict = FALSE, nonzero = FALSE)
}
\arguments{
  \item{i,j,k}{
    Integer vectors of equal length (or length 1), specifying the cells in the
    array which have non-zero entries.
  }
  \item{x}{
    Vector (numeric, integer, logical or complex)
    of the same length as \code{i}, \code{j} and \code{k},
    giving the values of the array entries that are not zero.
  }
  \item{dims}{
    Dimension of the array. An integer vector of length 3.
  }
  \item{dimnames}{
    Names for the three margins of the array. Either \code{NULL}
    or a list of three character vectors.
  }
  \item{strict}{
    Logical value specifying whether to enforce the rule that
    each entry in \code{i,j,k,x} refers to a different cell.
    If \code{strict=TRUE}, entries which refer to the same cell
    in the array will be reduced to a single entry by summing the
    \code{x} values. Default is \code{strict=FALSE}.
  }
  \item{nonzero}{
    Logical value specifying whether to remove any entries
    of \code{x} which equal zero.
  }
}
\details{
  An array \code{A} is three-dimensional if it is indexed by
  three integer indices, so that \code{A[i,j,k]} specifies an element of
  the array. The array is called sparse if only a small fraction of the
  entries are non-zero. A sparse array can be represented economically
  by listing only the entries which are non-zero.

  The \pkg{spatstat.sparse} package defines the class
  \code{sparse3Darray} of sparse three-dimensional arrays.
  These arrays can have numeric, integer, logical, or complex
  entries.

  The function \code{sparse3Darray} creates an object of class
  \code{"sparse3Darray"}. This object is essentially a list containing
  the vectors \code{i,j,k,x} and the arguments \code{dims,dimnames}.

  The arguments \code{i,j,k,x} should be vectors of equal length 
  identifying the cells in the array which have non-zero entries (indexed
  by \code{i,j,k}) and giving the values in these cells (given by \code{x}).

  The default behaviour of \code{sparse3Darray}
  is to accept the arguments \code{i,j,k,x} without modifying them.
  This would allow some entries of \code{x} to be equal to zero,
  and would allow a cell in the array to be referenced more than once
  in the indices \code{i,j,k}.

  If \code{nonzero=TRUE}, entries will be removed if the \code{x} value
  equals zero.

  If \code{strict=TRUE}, entries which refer to the same cell in the
  array will be combined into a single entry by summing the \code{x} values.
}
\value{
  An object of class \code{"sparse3Dvector"}.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{as.sparse3Darray}}
}
\examples{
  ## creation by specifying nonzero elements
  M <- sparse3Darray(i=1:3, j=c(3,1,2), k=4:2,
                     x=runif(3), dims=rep(4, 3))
  M
  ## duplicate entries
  Mn <- sparse3Darray(i=c(1,1,2), j=c(2,2,1), k=c(3,3,2),
                     x=runif(3), dims=rep(3, 3))
  ## cumulate entries in duplicate positions
  Ms <- sparse3Darray(i=c(1,1,2), j=c(2,2,1), k=c(3,3,2),
                      x=runif(3), dims=rep(3, 3), strict=TRUE)
}
\keyword{array}
\concept{sparse}
