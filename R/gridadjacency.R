#'
#'   gridadjacency.R
#'
#'   Adjacency matrix for points on a 2D integer grid.
#'
#'   $Revision: 1.1 $ $Date: 2021/03/15 05:59:26 $
#'

gridadjacencymatrix <- function(dims, across=TRUE, down=TRUE, diagonal=TRUE) {
  dims <- ensure2vector(dims)
  nr <- dims[1]
  nc <- dims[2]
  n <- prod(dims)
  serial <- matrix(1:n, nr, nc)
  m <- sparseMatrix(i=integer(0), j=integer(0), x=logical(0), dims=c(n,n))
  if(across) {
    #' join cells in adjacent columns (i, j) ~ (i, j+1) across each row
    allbutlastcol  <- as.vector(serial[   , -nc, drop=FALSE])
    allbutfirstcol <- as.vector(serial[   ,  -1, drop=FALSE])
    m[cbind(allbutfirstcol, allbutlastcol)] <- TRUE
  }
  if(down) {
    #' join cells in adjacent rows (i, j) ~ (i+1, j) down each column
    allbutlastrow  <- as.vector(serial[-nr,    , drop=FALSE])
    allbutfirstrow <- as.vector(serial[ -1,    , drop=FALSE])
    m[cbind(allbutfirstrow, allbutlastrow)] <- TRUE
  }
  if(diagonal) {
    #' join cells (i, j) ~ (i+1, j+1) 
    allexcbotleft  <-  as.vector(serial[-1,   -1, drop=FALSE])
    allexctopright <-  as.vector(serial[-nr, -nc, drop=FALSE])
    m[cbind(allexcbotleft, allexctopright)] <- TRUE
    #' join cells (i, j) ~ (i+1, j-1) 
    allexcbotright <-  as.vector(serial[-nr, -1, drop=FALSE])
    allexctopleft  <-  as.vector(serial[-1, -nc, drop=FALSE])
    m[cbind(allexcbotright, allexctopleft)] <- TRUE
  }
  m <- m | t(m)
  return(m)
}


