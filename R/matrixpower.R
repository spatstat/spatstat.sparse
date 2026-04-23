#'
#'       matrixpower.R
#'
#'   Copyright (c) Adrian Baddeley, Ege Rubak and Rolf Turner 2016-2022
#'   GNU Public Licence >= 2.0
#'
#'   $Revision: 1.9 $  $Date: 2024/06/09 00:01:09 $
#' 

matrixsqrt <- function(x, complexOK=TRUE) {
  ## matrix square root
  if(length(dim(x)) != 2)
    stop("x must be a matrix")
  if(!is.matrix(x))
    x <- as.matrix(x)
  check.nmatrix(x, mname="x") ## requires square matrix
  if(missing(complexOK) && is.complex(x)) complexOK <- TRUE
  if(!complexOK) stopifnot(is.numeric(x)) else
                 stopifnot(is.numeric(x) || is.complex(x))
  ## eigen decomposition
  e <- eigen(x)
  values <- e$values
  vectors <- e$vectors
  ## real or complex?
  realresult <- is.numeric(vectors) && is.numeric(values) && all(values >= 0)
  if(realresult) {
    ## result is a real matrix
    y <- vectors %*% diag(sqrt(values)) %*% t(vectors)
  } else {
    ## result is a complex matrix
    if(!complexOK)                                              
      stop("square root is complex, but complexOK=FALSE was specified", call.=FALSE)
    values <- as.complex(values) 
    y <- vectors %*% diag(sqrt(values)) %*% solve(vectors)
  }
  ## add dimnames
  if(!is.null(dn <- dimnames(x)))
    dimnames(y) <- rev(dn)
  return(y)
}

matrixinvsqrt <- function(x, complexOK=TRUE) {
  ## matrix inverse square root
  if(length(dim(x)) != 2)
    stop("x must be a matrix")
  if(!is.matrix(x))
    x <- as.matrix(x)
  check.nmatrix(x, mname="x") ## requires square matrix
  if(missing(complexOK) && is.complex(x)) complexOK <- TRUE
  if(!complexOK) stopifnot(is.numeric(x)) else
                 stopifnot(is.numeric(x) || is.complex(x))
  ## eigen decomposition
  e <- eigen(x)
  values <- e$values
  vectors <- e$vectors
  if(any(values == 0))
    stop("matrix is singular; cannot compute inverse square root", call.=FALSE)
  ## real or complex?
  realresult <- is.numeric(vectors) && is.numeric(values) && all(values >= 0)
  if(realresult) {
    ## result is a real matrix
    y <- vectors %*% diag(1/sqrt(values)) %*% t(vectors)
  } else {
    ## result is a complex matrix
    if(!complexOK)
      stop("inverse square root is complex, but complexOK=FALSE was specified", call.=FALSE)
    values <- as.complex(values) 
    y <- vectors %*% diag(1/sqrt(values)) %*% solve(vectors)
  }
  ## add dimnames
  if(!is.null(dn <- dimnames(x)))
    dimnames(y) <- rev(dn)
  return(y)
}

matrixpower <- function(x, power, complexOK=TRUE) {
  check.1.real(power)
  if(length(dim(x)) != 2)
    stop("x must be a matrix")
  if(!is.matrix(x))
    x <- as.matrix(x)
  check.nmatrix(x, mname="x") ## requires square matrix
  if(power == 0) {
    ## power = 0 yields identity matrix even if x is singular
    y <- diag(nrow(x))
  } else {
    ## nonzero power
    if(missing(complexOK) && is.complex(x)) complexOK <- TRUE
    if(!complexOK) stopifnot(is.numeric(x)) else
                   stopifnot(is.numeric(x) || is.complex(x))
    ## eigen decomposition
    e <- eigen(x)
    values <- e$values
    vectors <- e$vectors
    if(any(values == 0) && power < 0)
      stop("matrix is singular; cannot compute negative power", call.=FALSE)
    ## real or complex?
    realresult <- is.numeric(vectors) && is.numeric(values) && all(values >= 0)
    if(realresult) {
      ## result is a real matrix
      y <- vectors %*% diag(values^power) %*% t(vectors)
    } else {
      ## result is a complex matrix
      if(!complexOK)
        stop("result is complex, but complexOK=FALSE was specified", call.=FALSE)
      values <- as.complex(values)
      y <- vectors %*% diag(values^power) %*% solve(vectors)
    }
  }
  ## add dimnames
  if(!is.null(dn <- dimnames(x)))
    dimnames(y) <- rev(dn)
  return(y)
}

symmatrix <- function(x, from=c("lower", "upper"), diag=TRUE) {
  from <- match.arg(from)
  x <- as.vector(x)
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
  M <- matrix(, m, m)
  A <- switch(from,
              lower = lower.tri(M, diag=diag),
              upper = upper.tri(M, diag=diag))
  M[  A ] <- x
  M[ !A ] <- t(M)[ !A ]
  return(M)
}

