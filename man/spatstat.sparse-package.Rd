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
  These arrays can have numeric, integer, factor, logical, or complex
  entries. 
  
  The package supports:
  \itemize{
    \item creation of sparse arrays, conversion to/from other data types
    \item extracting of entries
    \item assigning new values to entries
    \item permuting and combining arrays
    \item arithmetic and logical operations
    \item tensor operations (generalising matrix multiplication)
    \item printing of sparse arrays.
  }

  The \pkg{spatstat.sparse} package uses the \pkg{Matrix} package to
  handle slices of three-dimensional arrays which are two-dimensional
  (sparse matrices) or one-dimensional (sparse vectors).

  Most functions are currently in \code{spatstat.sparse-internal}.
  Here is the full list:

  \itemize{
    \item \code{sparse3Darray}
    \item \code{as.sparse3Darray}
    \item \code{dim.sparse3Darray}
    \item \code{"dim<-.sparse3Darray"}
    \item \code{dimnames.sparse3Darray}
    \item \code{"dimnames<-.sparse3Darray"}
    \item \code{print.sparse3Darray}
    \item \code{aperm.sparse3Darray}
    \item \code{as.array.sparse3Darray}
    \item \code{"[.sparse3Darray"}
    \item \code{"[<-.sparse3Darray"}
    \item \code{anyNA.sparse3Darray}
    \item \code{RelevantZero}
    \item \code{RelevantNA}
    \item \code{RelevantEmpty}
    \item \code{isRelevantZero}
    \item \code{unionOfSparseIndices}
    \item \code{Complex.sparse3Darray}
    \item \code{Math.sparse3Darray}
    \item \code{Ops.sparse3Darray}
    \item \code{Summary.sparse3Darray}
    \item \code{inside3Darray}
    \item \code{SparseEntries}
    \item \code{SparseIndices}
    \item \code{EntriesToSparse}
    \item \code{mapSparseEntries}
    \item \code{applySparseEntries}
    \item \code{sumsymouterSparse}
    \item \code{tensorSparse}
    \item \code{marginSumsSparse}
    \item \code{rbindCompatibleDataFrames}
    \item \code{bind.sparse3Darray}
    \item \code{check.anySparseVector}
    \item \code{evalSparse3Dentrywise}
    \item \code{expandSparse}
    \item \code{sparseVectorCumul}
    \item \code{tensor1x1}
  }
  Indexing utilities:
  \itemize{
    \item \code{fullIndexSequence}
    \item \code{grokIndexVector}
    \item \code{replacementIndex}
    \item \code{representativeRows}
    \item \code{strictIndexSequence}
  }
}
\section{Matrix Utilities}{
Existing help files:
  \code{\link{sumouter}}, \code{\link{quadform}},
  \code{\link{bilinearform}}.

  Functions to document:
  \itemize{
    \item \code{sumsymouter}
    \item \code{checksolve}
    \item \code{check.mat.mul}
  }
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
