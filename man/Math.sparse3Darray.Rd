\name{Math.sparse3Darray}
\alias{Math.sparse3Darray}
\alias{Ops.sparse3Darray}
\alias{Complex.sparse3Darray}
\alias{Summary.sparse3Darray}
\title{S3 Group Generic Methods for Sparse Three-Dimensional Arrays}
\description{
  Group generic methods which make it possible to
  apply the familiar mathematical operators and functions
  to sparse three-dimensional arrays (objects of class
  \code{"sparse3Darray"}).
  See Details for a list of implemented functions.
}
\usage{
## S3 methods for group generics have prototypes:
\special{Math(x, \dots)}
\special{Ops(e1, e2)}
\special{Complex(z)}
\special{Summary(\dots, na.rm=FALSE)}
%NAMESPACE S3method("Math", "sparse3Darray")
%NAMESPACE S3method("Ops", "sparse3Darray")
%NAMESPACE S3method("Complex", "sparse3Darray")
%NAMESPACE S3method("Summary", "sparse3Darray")
}
\arguments{
  \item{x, z, e1, e2}{
    Sparse three-dimensional arrays (objects of class
    \code{"sparse3Darray"}). Alternatively \code{e1} or \code{e2} can be
    a single scalar, vector, sparse vector, matrix or sparse matrix.
  }
  \item{\dots}{further arguments passed to methods.}
  \item{na.rm}{
    Logical value specifying whether missing values
    should be removed. 
  }
}
\details{
  These group generics make it possible to
  perform element-wise arithmetic and logical operations
  with sparse three-dimensional arrays, or apply mathematical functions
  element-wise, or compute standard summaries such as the mean and maximum.
  
  Below is a list of mathematical functions and operators which are
  defined for sparse 3D arrays. 

  \enumerate{
    \item Group \code{"Math"}:
    \itemize{
      \item
      \code{abs}, \code{sign}, \code{sqrt},\cr
      \code{floor}, \code{ceiling}, \code{trunc},\cr
      \code{round}, \code{signif}

      \item
      \code{exp}, \code{log},  \code{expm1}, \code{log1p},\cr
      \code{cos}, \code{sin}, \code{tan},\cr
      \code{cospi}, \code{sinpi}, \code{tanpi},\cr
      \code{acos}, \code{asin}, \code{atan}

      \code{cosh}, \code{sinh}, \code{tanh},\cr
      \code{acosh}, \code{asinh}, \code{atanh}

      \item
      \code{lgamma}, \code{gamma}, \code{digamma}, \code{trigamma}
      \item \code{cumsum}, \code{cumprod}, \code{cummax}, \code{cummin}
    }

    \item Group \code{"Ops"}:
    \itemize{
      \item
      \code{"+"}, \code{"-"}, \code{"*"}, \code{"/"},
      \code{"^"}, \code{"\%\%"}, \code{"\%/\%"}

      \item \code{"&"}, \code{"|"}, \code{"!"}

      \item \code{"=="}, \code{"!="},
      \code{"<"}, \code{"<="}, \code{">="}, \code{">"}
    }

    \item Group \code{"Summary"}:
    \itemize{
      \item \code{all}, \code{any}
      \item \code{sum}, \code{prod}
      \item \code{min}, \code{max}
      \item \code{range}
    }

    \item Group \code{"Complex"}:
    \itemize{
      \item \code{Arg}, \code{Conj}, \code{Im}, \code{Mod}, \code{Re}
    }
  }
}
\seealso{
  \code{\link{sparse3Darray}},
  \code{\link{tensorSparse}}
}
\examples{
  M <- sparse3Darray(i=1:4, j=sample(1:4, replace=TRUE),
                     k=c(1,2,1,2), x=1:4, dims=c(5,5,2))
  negM <- -M
  twoM <- M + M
  Mplus <- M + 1 ## not sparse!
  posM <- (M > 0)
  range(M)
  sinM <- sin(M)
  cosM <- cos(M) ## not sparse!
  expM1 <- expm1(M)
}
\value{
  The result of group \code{"Math"} functions is another 
  three-dimensional array of the same dimensions as \code{x},
  which is sparse if the function maps 0 to 0, and otherwise is a
  full three-dimensional array.

  The result of group \code{"Ops"} operators is 
  another three-dimensional array of the same dimensions as
  \code{e1} and \code{e2}, which is sparse if both \code{e1} and
  \code{e2} are sparse.

  The result of group \code{"Complex"} functions is 
  another sparse three-dimensional array of the same dimensions as
  \code{z}.

  The result of group \code{"Summary"} functions is
  a logical value or a numeric value or a numeric vector of length 2.
}
\author{
  \spatstatAuthors.
}
\keyword{spatial}
\keyword{methods}
