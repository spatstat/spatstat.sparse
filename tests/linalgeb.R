#' Header for spatstat.sparse/tests/*R
#'

require(spatstat.sparse)
ALWAYS <- FULLTEST <- TRUE
##
##    tests/linalgeb.R
##
## checks validity of linear algebra code
##
##  $Revision: 1.5 $ $Date: 2020/01/05 02:34:17 $
##

local({
  p <- 3
  n <- 4
  k <- 2
  
  x <- matrix(1:(n*p), n, p)
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
  
  x <- array(as.numeric(1:(p * n * n)), dim=c(p, n, n))
  w <- matrix(1:(n*n), n, n)
  y <- matrix(numeric(p * p), p, p)
  for(i in 1:n)
    for(j in (1:n)[-i])
      y <- y + w[i,j] * outer(x[,i,j], x[,j,i])
  z <- sumsymouter(x, w)
  if(!identical(y,z))
    stop("sumsymouter gives incorrect result")

  #' power of complex matrix
  M <- diag(c(1,-1))
  V <- matrixsqrt(M)
  V <- matrixinvsqrt(M)
  V <- matrixpower(M, 1/2)

  #' infrastructure
  A <- matrix(1:12, 3, 4)
  B <- matrix(1:8, 4, 2)
  check.mat.mul(A, B)
  check.mat.mul(A, B[,1])
  check.mat.mul(A, A, fatal=FALSE)
})
