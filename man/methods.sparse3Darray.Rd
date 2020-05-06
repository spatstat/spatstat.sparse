\name{methods.sparse3Darray}
\Rdversion{1.1}
\alias{methods.sparse3Darray} %DoNotExport
\alias{anyNA.sparse3Darray}
\alias{dim.sparse3Darray}
\alias{dim<-.sparse3Darray}
\alias{dimnames.sparse3Darray}
\alias{dimnames<-.sparse3Darray}
\alias{print.sparse3Darray}
\title{
  Methods for Sparse Three-Dimensional Arrays
}
\description{
  Methods for the class \code{"sparse3Darray"}
  of sparse three-dimensional arrays.
}
\usage{
 \method{anyNA}{sparse3Darray}(x, recursive = FALSE)
 \method{dim}{sparse3Darray}(x)
 \method{dim}{sparse3Darray}(x) <- value
 \method{dimnames}{sparse3Darray}(x)
 \method{dimnames}{sparse3Darray}(x) <- value
 \method{print}{sparse3Darray}(x, \dots)
}
\arguments{
  \item{x}{
    A sparse three-dimensional array (object of class \code{"sparse3Darray"}).
  }
  \item{value}{
    Replacement value (see Details).
  }
  \item{recursive,\dots}{
    Ignored.
  }
}
\details{
  These are methods for the generics
  \code{\link[base]{anyNA}},
  \code{\link[base]{dim}},
  \code{\link[base]{dim<-}},
  \code{\link[base]{dimnames}},
  \code{\link[base]{dimnames<-}}
  and \code{\link[base]{print}}
  for the class \code{"sparse#Darray"} of sparse three-dimensional
  arrays.

  For \code{dimnames(x) <- value}, the \code{value} should either be
  \code{NULL}, or a list of length 3 containing character vectors giving
  the names of the margins.
  
  For \code{dim(x) <- value}, the \code{value} should be an integer
  vector of length 3 giving the new dimensions of the array. Note that this
  operation does not change the array positions of the non-zero entries
  (unlike \code{dim(x) <- value} for a full array). An error occurs if
  some of the non-zero entries would lie outside the new extent of the array.
}
\value{
  \code{anyNA} returns a single logical value.

  \code{dim} returns an integer vector of length 3.

  \code{dimnames} returns \code{NULL}, or a list of length 3
  whose entries are character vectors.

  \code{dim<-} and \code{dimnames<-} return a sparse 3D array.

  \code{print} returns \code{NULL}, invisibly.
}
\author{
  \spatstatAuthors.
}
\seealso{
  \code{\link{sparse3Darray}}
}
\examples{
 M <- sparse3Darray(i=1:4, j=sample(1:4, replace=TRUE),
                    k=c(1,2,1,2), x=1:4, dims=c(5,5,2))

 anyNA(M)
 dim(M)
 dimnames(M)
 dimnames(M) <- list(letters[1:5], LETTERS[1:5], c("Yes", "No"))
 print(M)
}
\keyword{array}
\keyword{manip}
\concept{sparse}
