\name{Extract.sparse3Darray}
\alias{[.sparse3Darray}
\alias{[<-.sparse3Darray}
\title{Extract or Replace Entries in a Sparse Array}
\description{
  Extract or replace entries in a sparse three-dimensional array.
}
\usage{
  \method{[}{sparse3Darray}(x, i, j, k, drop=TRUE, \dots)
  \method{[}{sparse3Darray}(x, i, j, k, \dots) <- value
}
\arguments{
  \item{x}{
    Sparse three-dimensional array
    (object of class \code{"sparse3Darray"}).
  }
  \item{i,j,k}{
   Subset indices for each dimension of the array. See Details.
  }
  \item{value}{
    Replacement value for the subset. 
  }
  \item{drop}{
    Logical value indicating whether to return a lower-dimensional
    object (matrix or vector) when appropriate.
  }
  \item{\dots}{
    Ignored. This argument is required for compatibility
    with the generic function.
  }
}
\value{
  \code{[.sparse3Darray} returns either
  a sparse three-dimensional array (class \code{"sparse3Darray"}),
  a sparse matrix (class \code{sparseMatrix} in the \pkg{Matrix} package),
  a sparse vector (class \code{sparseVector} in the \pkg{Matrix} package),
  or in some cases a full array, matrix or vector.

  \code{[<-.sparse3Darray} returns another sparse three-dimensional
  array.
}
\details{
  These functions are defined for a sparse three-dimensional array \code{x}.
  They extract a designated subset of the array,
  or replace the values in the designated subset.

  The function \code{[.sparse3Darray} is a method for the generic
  subset extraction operator \code{\link{[}}.
  The function \code{[<-.sparse3Darray} is a method for the generic
  subset replacement operator \code{\link{[<-}}.

  These methods use the same indexing rules
  as the subset operator for full arrays:
  \itemize{
    \item If \code{i}, \code{j} and \code{k} are integer vectors,
    the subset is the Cartesian product (i.e. all cells in the array
    identified by an entry of \code{i}, an entry of \code{j} and
    an entry of \code{k}).
    \item Some or all of the arguments \code{i}, \code{j} and \code{k}
    may be missing from the call; a missing index argument is
    interpreted as meaning that all possible values of that index are
    allowed.
   \item Arguments \code{i}, \code{j} and \code{k} may be logical
     vectors (with the value \code{TRUE} assigned to entries that should
     be included).
   \item Arguments \code{i}, \code{j} and \code{k} may be
     character vectors with entries matching the corresponding
     \code{dimnames}.
   \item Argument \code{i} may be an integer matrix with 3 columns
      (and the arguments \code{j,k} should be absent).
      Each row of the matrix contains the indices
      of one cell in the array.
   }

   If the designated subset lies within the array bounds, then
   the result of \code{[} will be a sparse three-dimensional array,
   sparse matrix or sparse vector. If \code{drop=FALSE} the result
   will always be three-dimensional; if \code{drop=TRUE} (the default)
   the result will be reduced to two or one dimensions when appropriate.

   If the designated subset \emph{does not} lie within the array bounds, then
   the result of \code{[} will be a full three-dimensional array,
   matrix or vector containing \code{NA} values at the positions that
   were outside the array bounds.

   The result of \code{[<-} is always a sparse three-dimensional array.
   If the designated subset did not lie within the array bounds of
   \code{x}, then the array bounds will be extended (with a warning message).
}
\seealso{
  \code{\link{sparse3Darray}},
  \code{\link{methods.sparse3Darray}}.
}
\examples{
  M <- sparse3Darray(i=1:4, j=sample(1:4, replace=TRUE),
                     k=c(1,2,1,2), x=1:4, dims=c(5,5,2))
  dimnames(M) <- list(letters[1:5], LETTERS[1:5], c("yes", "no"))
  M[ 3:4, , ]
  M[ 3:4, 2:4, ]
  M[ 4:3, 4:2, 1:2]
  M[, 3, ]
}
\author{
  \spatstatAuthors.
}
\keyword{array}
\keyword{manip}
\concept{sparse}
