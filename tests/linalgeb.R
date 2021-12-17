#' Header for spatstat.sparse/tests/*R
#'

require(spatstat.sparse)
ALWAYS <- FULLTEST <- TRUE
##
##    tests/linalgeb.R
##
## checks validity of linear algebra code
##
##  $Revision: 1.12 $ $Date: 2021/12/17 00:32:32 $
##

local({
  p <- 3
  n <- 4
  k <- 2

  ## correctness of 'quadform'
  x <- matrix(1:(n*p), n, p)
  v <- matrix(runif(p^2), p, p)
  zW <- zU <- numeric(n)
  for(i in 1:n) {
    xi <- x[i,,drop=FALSE]
    zW[i] <- xi %*% v %*% t(xi)
    zU[i] <- xi %*% t(xi)
  }
  if(!all(zU == quadform(x)))
    stop("quadform gives incorrect result in Unweighted case")
  if(!all(zW == quadform(x,v)))
    stop("quadform gives incorrect result in Weighted case")
    
  ## correctness of 'sumouter'
  w <- runif(n)
  y <- matrix(1:(2*n), n, k)
  zUS <- zWS <- matrix(0, p, p)
  zUA <- zWA <- matrix(0, p, k)
  for(i in 1:n) {
    zUS <- zUS +        outer(x[i,],x[i,])
    zWS <- zWS + w[i] * outer(x[i,],x[i,])
    zUA <- zUA +        outer(x[i,],y[i,])
    zWA <- zWA + w[i] * outer(x[i,],y[i,])
  }
  if(!identical(zUS, sumouter(x)))
    stop("sumouter gives incorrect result in Unweighted Symmetric case")
  if(!identical(zWS, sumouter(x,w)))
    stop("sumouter gives incorrect result in Weighted Symmetric case")
  if(!identical(zUA, sumouter(x, y=y)))
    stop("sumouter gives incorrect result in Unweighted Asymmetric case")
  if(!identical(zWA, sumouter(x, w, y)))
    stop("sumouter gives incorrect result in Weighted Asymmetric case")

  #' complex quadratic forms - execute only
  dimnames(x) <- list(letters[1:nrow(x)], LETTERS[1:ncol(x)])
  a <- quadform(x + 1i)
  b <- quadform(x + 1i, v+1i)
  d <- quadform(x     , v+1i)
  a <- sumouter(x + 1i)
  b <- sumouter(x + 1i, w + 1i)
  d <- sumouter(x + 1i, w + 1i, x - 1i)

  #' NA values
  xna <- x; xna[1,1] <- NA
  wna <- w; w[2] <- NA
  vna <- v; v[1,2] <- NA
  o <- quadform(xna)
  o <- quadform(xna, vna)
  o <- sumouter(xna)
  o <- sumouter(xna, w)
  o <- sumouter(xna, wna)
  
  #' sumsymouter
  x <- array(as.numeric(1:(p * n * n)), dim=c(p, n, n))
  w <- matrix(1:(n*n), n, n)
  y <- matrix(numeric(p * p), p, p)
  #' check correctness
  for(i in 1:n)
    for(j in (1:n)[-i])
      y <- y + w[i,j] * outer(x[,i,j], x[,j,i])
  z <- sumsymouter(x, w)
  if(!identical(y,z))
    stop("sumsymouter gives incorrect result")
  #' cover code blocks
  o <- sumsymouter(x, distinct=FALSE)
  a <- sumsymouter(x + 1i)
  b <- sumsymouter(x + 1i, w + 1i)
  if(require(Matrix)) o <- sumsymouter(x, as(w, "sparseMatrix"))
  
  #' power of complex matrix
  M <- diag(c(4,-4))
  dimnames(M) <- list(letters[1:2], letters[1:2])
  V <- matrixsqrt(M)
  V <- matrixinvsqrt(M)
  V <- matrixpower(M, 1/2)
  U <- matrixsqrt(abs(M), complexOK=FALSE)

  #' infrastructure
  A <- matrix(1:12, 3, 4)
  B <- matrix(1:8, 4, 2)
  check.mat.mul(A, B)
  check.mat.mul(A, B[,1])
  check.mat.mul(A, A, fatal=FALSE)
  D <- diag(c(1,4,9))
  checksolve(D)
  D[1,1] <- 0
  checksolve(D, "silent")
})
