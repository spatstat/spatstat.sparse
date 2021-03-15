\name{gridadjacencymatrix}
\alias{gridadjacencymatrix}
\title{
  Create Adjacency Matrix for Spatial Grid
}
\description{
  Given the dimensions of a rectangular grid of points, this command
  creates the adjacency matrix for the corresponding neighbourhood
  graph, whose vertices are the grid points, and whose edges are the
  joins between neighbouring grid points.
}
\usage{
gridadjacencymatrix(dims, across = TRUE, down = TRUE, diagonal=TRUE)
}
\arguments{
  \item{dims}{
    Grid dimensions. An integer, or a vector of two integers.
    First entry specifies the number of points in the \eqn{y} direction.
  }
  \item{across}{
    Logical value equal to \code{TRUE} if horizontal neighbours should
    be joined.
  }
  \item{down}{
    Logical value equal to \code{TRUE} if vertical neighbours should
    be joined.
  }
  \item{diagonal}{
    Logical value equal to \code{TRUE} if diagonal neighbours should
    be joined.
  }
}
\details{
  If \code{N = prod(dims)} is the total number of grid points,
  then the result is an \code{N * N} sparse matrix with logical entries
  equal to \code{TRUE} if the corresponding grid points are joined.
}
\value{
  A sparse matrix. 
}
\author{
  Adrian Baddeley.
}
\examples{
  gridadjacencymatrix(c(2,3))  
}
\keyword{datagen}
