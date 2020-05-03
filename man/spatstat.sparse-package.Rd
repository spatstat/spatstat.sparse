\name{spatstat.sparse-package}
\alias{spatstat.sparse-package} 
\alias{spatstat.sparse} 
\docType{package}
\title{The spatstat.sparse Package}
\description{
  The \pkg{spatstat.sparse} package defines three-dimensional
  sparse arrays, and supports standard operations on them.
  It also provides some utility functions for matrix calculations
  that are common in statistics, such as quadratic forms.
}
\details{
  The \pkg{spatstat.sparse} package
  \itemize{
    \item defines a class of sparse three-dimensional arrays
    and supports standard operations on them;
    \item provides utility functions for matrix computations
    such as quadratic forms.
  }
  The main purpose of \pkg{spatstat.sparse} is to 
  define a class of sparse three-dimensional arrays.

  An array \code{A} is three-dimensional if it is indexed by
  three integer indices, so that \code{A[i,j,k]} specifies an element of
  the array. The array is called sparse if only a small fraction of the
  entries are non-zero. A sparse array can be represented economically
  by listing only the entries which are non-zero.

  The \pkg{spatstat.sparse} package defines the class
  \code{sparse3Darray} of sparse three-dimensional arrays.
  These arrays can have numeric, integer, factor, logical, or complex
  entries. 
  
  The package supports:
  \itemize{
    \item extracting of entries
    \item assigning new values to entries
    \item permuting and combining arrays
    \item arithmetic and logical operations
    \item tensor operations (generalising matrix multiplication)
    \item conversion to other data types
    \item printing of sparse arrays.
  }

  The \pkg{spatstat.sparse} package uses the \pkg{Matrix} package to
  handle slices of three-dimensional arrays which are two-dimensional
  (sparse matrices) or one-dimensional (sparse vectors).

  MOST FUNCTIONS ARE CURRENTLY IN spatstat.sparse-internal.

  Also mention \code{\link{sumouter}}, \code{\link{quadform}},
  \code{\link{bilinearform}}.
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
