\name{runSparseMarkovChain}
\alias{runSparseMarkovChain}
\title{
  Run a Markov Chain with a Sparse Transition Matrix
}
\description{
  Generate a simulated realisation of a discrete-time finite-state
  Markov chain, using a sparse matrix representation of the transition
  probability matrix.
}
\usage{
  runSparseMarkovChain(P, x0, nsteps, \dots,
        result = c("history", "last"), check = TRUE,
        method = c("C", "interpreted"))
}
\arguments{
  \item{P}{
    Transition matrix. A sparse matrix supported by the
    \pkg{Matrix} package, or a \code{matrix} in the base \R package.
  }
  \item{x0}{
    Integer between \code{1} and \code{nrow(P)}
    specifying the initial state. Alternatively a vector of integers.
  }
  \item{nsteps}{
    Integer. Number of steps of the Markov chain.
  }
  \item{\dots}{
    Ignored.
  }
  \item{result}{
    Character string (partially matched) indicating whether to return
    the entire trajectory of the chain (\code{result="history"}, the
    default) or just the final state of the chain (\code{result="last"}).
  }
  \item{check}{
    Logical value specifying whether to check the validity of \code{P}
    and \code{x0}.
  }
  \item{method}{
    Character string (partially matched) indicating whether to use
    the \proglang{C} code implementation or an interpreted \R
    implementation. For testing purposes only.
  }
}
\details{
  If \code{x0} is a single integer, the Markov chain
  with transition probability matrix \code{P} will be
  run for \code{nsteps} steps starting from initial state \code{x0}.
  The result will be a vector of length \code{nsteps+1} giving the
  history of the chain after \code{0, 1, \dots, nsteps} steps
  (if \code{result="history"}, the default) or a
  single integer giving the final state of the chain
  (if \code{result="last"}).

  If \code{x0} is a vector of integers, each entry of \code{x0}
  will be taken as the starting state for an independent copy of the
  Markov chain. Each copy will be run for \code{nsteps} steps.
  The result will be a matrix with \code{n=length(x0)} rows
  and \code{nsteps+1} columns giving the history of each chain
  (if \code{result="history"}, the default) or a
  vector of \code{n} integers giving the final state of each chain
  (if \code{result="last"}).

  The matrix \code{P} should ideally be a sparse matrix
  in row major form (class \code{"dgRMatrix"}) but will be converted
  to this form.
}
\value{
  Integer vector or matrix.
}
\author{
  \adrian
}
\examples{
  M <- matrix(c(1, 0, 1, 0, 0,
                0, 1, 1, 0, 0,
                0, 1, 1, 1, 0,
                1, 1, 1, 1, 1,
                1, 1, 1, 1, 1),
              5, 5) 
  P <- M/rowSums(M)
  runSparseMarkovChain(P, 2, 20)

  ## random walk on { 1, ..., 1000 }
  G <- diag(1000)
  G[abs(row(G) - col(G)) == 1] <- 1
  H <- G/rowSums(G)
  plot(runSparseMarkovChain(H, 500, 1e4), type="l")
}
\keyword{datagen}
