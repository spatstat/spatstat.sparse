\name{spatstat.sparse-package}
\alias{spatstat.sparse-package} 
\alias{spatstat.sparse} 
\docType{package}
\title{The spatstat.sparse Package}
\description{
  The \pkg{spatstat.sparse} package defines three-dimensional
  sparse arrays, and supports standard operations on them.
  It also provides some utility functions for matrix calculations
  such as quadratic forms.
}
\details{
  The \pkg{spatstat.sparse} package
  \itemize{
    \item defines a class of sparse three-dimensional arrays
    and supports standard operations on them
    (see Section \emph{Sparse 3D Arrays}).
    \item provides utility functions for matrix computations
    that are common in statistics,
    such as quadratic forms (see Section \emph{Matrix Utilities}).
  }
  The code in \pkg{spatstat.sparse} was originally written
  for internal use within the \pkg{spatstat} package, but has now been
  removed and organised into a separate, stand-alone package
  which can be used for other purposes.
}
\section{Sparse 3D Arrays}{
  The main purpose of \pkg{spatstat.sparse} is to 
  define a class of sparse three-dimensional arrays.

  An array \code{A} is three-dimensional if it is indexed by
  three integer indices, so that \code{A[i,j,k]} specifies an element of
  the array. The array is called sparse if only a small fraction of the
  entries are non-zero. A sparse array can be represented economically
  by listing only the entries which are non-zero.

  The \pkg{spatstat.sparse} package defines the class
  \code{sparse3Darray} of sparse three-dimensional arrays.
  These arrays can have numeric, integer, logical, or complex
  entries. 
  
  The package supports:
  \itemize{
    \item creation of sparse arrays from raw data
    \item conversion to/from other data types
    \item array indexing, extraction of entries, assignment of new values
    \item arithmetic and logical operations
    \item tensor operations (generalising matrix multiplication)
    \item permutation of array dimensions
    \item binding of several arrays into a single array
    \item printing of sparse arrays.
  }

  The \pkg{spatstat.sparse} package uses the \pkg{Matrix} package to
  handle slices of three-dimensional arrays which are two-dimensional
  (sparse matrices) or one-dimensional (sparse vectors).

  The main functions are:
  \tabular{ll}{
    \code{\link{sparse3Darray}} \tab Create a sparse 3D array \cr
    \code{\link{as.sparse3Darray}} \tab Convert other data to a sparse
    3D array \cr
    \code{\link{[.sparse3Darray}} \tab Subset operator \cr
    \code{\link{aperm.sparse3Darray}} \tab Permute a sparse array \cr
    \code{\link{Ops.sparse3Darray}} \tab arithmetic and logical operators \cr
    \code{\link{Complex.sparse3Darray}} \tab complex operators \cr
    \code{\link{Math.sparse3Darray}}\tab standard mathematical functions \cr
    \code{\link{Summary.sparse3Darray}} \tab mean, maximum etc \cr
    \code{\link{tensorSparse}} \tab Tensor product \cr
    \code{\link{as.array.sparse3Darray}} \tab Convert sparse array to full array \cr
  }    

  The class \code{"sparse3Darray"} has
  methods for \code{anyNA}, \code{dim}, \code{dim<-},
  \code{dimnames}, \code{dimnames<-} and \code{print},
  documented in \code{\link{methods.sparse3Darray}}.

  For other undocumented functions,
  see \code{\link{spatstat.sparse-internal}}.
}
\section{Matrix Utilities}{
  The package also includes some utilities for matrix calculations:
  \tabular{ll}{
    \code{\link{sumouter}} \tab sum of outer products of rows of a matrix \cr
    \code{\link{quadform}} \tab quadratic form involving rows of a matrix \cr
    \code{\link{bilinearform}} \tab bilinear form involving rows of a matrix \cr
    \code{\link{matrixsqrt}} \tab square root of a matrix \cr
    \code{\link{matrixpower}} \tab powers of a matrix \cr
  }
  % Functions to document: sumsymouter, checksolve, check.mat.mul
}
\section{Licence}{
  This library and its documentation are usable under the terms of the
  \dQuote{GNU General Public License},
  a copy of which is distributed with \R.
}
\author{
  \spatstatAuthors.
}
\keyword{array}
\keyword{algebra}
\keyword{package}
