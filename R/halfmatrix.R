#'
#'  halfmatrix.R
#'
#'  Extract or reinstate lower/upper triangle of a matrix
#'
#'   Copyright (c) Adrian Baddeley, Ege Rubak and Rolf Turner 2026
#'   GNU Public Licence >= 2.0
#' 
#'  $Revision: 1.2 $ $Date: 2026/04/27 04:18:57 $

symmatrix <- function(x, from=c("lower", "upper"), diag=TRUE, na.value=NA) {
  from <- match.arg(from)
  ## infer dimensions
  n <- length(x)
  if(diag) {
    m <- (sqrt(8*n+1)-1)/2
    E <- m^2 + m - 2*n
  } else {
    m <- (sqrt(8*n+1)+1)/2
    E <- m^2 - m - 2*n
  }
  if(E != 0)
    stop("x has the wrong length for a triangular subset", call.=FALSE)
  ## go
  if(inherits(x, "sparseVector")) {
    k <- x@i
    switch(from,
           lower = {
             colHeights <- if(diag) m:1 else (m-1):0
             starts <- c(0, cumsum(colHeights))
             j <- findInterval(k, starts+0.5, all.inside=TRUE)
             i <- j - diag + (k - starts[j])
             M <- sparseMatrix(i=i, j=j, x=x@x, symmetric=TRUE,
                               dims=c(m,m))
           },
           upper = {
             rowWidths <- if(diag) m:1 else (m-1):0
             starts <- c(0, cumsum(rowWidths))
             i <- findInterval(k, starts+0.5, all.inside=TRUE)
             j <- i - diag + (k - starts[j])
             M <- sparseMatrix(i=i, j=j, x=x@x, symmetric=TRUE,
                               dims=c(m,m))
           })
    if(!diag) diag(M) <- na.value
  } else {
    x <- as.vector(x)
    M <- matrix(na.value, m, m)
    A <- switch(from,
                lower = lower.tri(M, diag=diag),
                upper = upper.tri(M, diag=diag))
    M[  A ] <- x
    M[ !A ] <- t(M)[ !A ]
  }
  return(M)
}

