\name{sumsymouter}
\alias{sumsymouter}
\title{Compute Quadratic Forms in 3D Arrays}
\description{
  Calculates certain quadratic forms of 3D arrays.
}
\usage{
  sumsymouter(x, w = NULL, distinct = TRUE)
}
\arguments{
  \item{x}{
    Three-dimensional array (class \code{array})
    or sparse three-dimensional array (class \code{sparse3Darray}).
  }
  \item{w}{Matrix or sparse matrix.}
  \item{distinct}{
    Logical value indicating whether to exclude diagonal contributions.
    See Details.
  }
}
\details{
  Given a three-dimensional matrix or sparse matrix \eqn{x}
  and a matrix or sparse matrix \eqn{w},
  this function computes the matrix \eqn{M} with entries
  \deqn{
    M_{i,j} = \sum_k \sum_m w_{k,m} x_{i, k, m} x_{j, m, k}
  }{
    M[i,j] = sum[k] sum[m] w[k,m] * x[i,k,m] * x[j, m, k]
  },
  the weighted sum of symmetric outer products of \eqn{x}
  with weights \eqn{w}.

  If \code{distinct=TRUE} (the default), terms with \eqn{k=m} are
  omitted from the sum. If \code{distinct=FALSE} they are included.

  The array \code{x} should have dimensions
  \eqn{p \times n \times n}{p * n * n}
  and the matrix \code{w} should have dimensions \eqn{n \times n}{n * n}
  if it is given. 
  The result \eqn{M} is a \eqn{p \times p}{p * p} square matrix.

  If \code{w} is missing or \code{NULL}, all weights are taken to be 1.

  Entries in \code{x} and \code{w} may be integer, numeric, complex
  or logical values.
  
  The result is a symmetric matrix if \code{w} is a symmetric matrix,
  or if \code{w} is missing or \code{NULL}. Otherwise it is not
  guaranteed to be symmetric.
}
\value{
  A square matrix with numeric or complex values.
}
\author{
  \adrian.
}
\seealso{
   \code{\link{sumouter}}
}
\examples{
   d <- c(2, 3, 3)
   n <- prod(d)
   v <- ifelse(runif(n) < 0.6, 0, round(runif(n), 1))
   x <- array(v, dim=d)
   w <- outer(1:3, 1:3)
   sumsymouter(x, w)
}
\keyword{array}
\keyword{algebra}
\concept{sparse}
