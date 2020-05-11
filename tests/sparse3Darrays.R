#' Header for spatstat.sparse/tests/*R
#'

require(spatstat.sparse)
ALWAYS <- FULLTEST <- TRUE
#'    tests/sparse3Darrays.R
#'  Basic tests of code in sparse3Darray.R and sparsecommon.R
#'  $Revision: 1.27 $ $Date: 2020/05/11 01:42:58 $

if(!exists("ALWAYS")) ALWAYS <- TRUE
if(!exists("FULLTEST")) FULLTEST <- ALWAYS

if(ALWAYS) { # fundamental, C code
local({
  #' forming arrays

  #' creation by specifying nonzero elements
  M <- sparse3Darray(i=1:3, j=c(3,1,2), k=4:2,
                     x=runif(3), dims=rep(4, 3))
  #' duplicate entries
  Mn <- sparse3Darray(i=c(1,1,2), j=c(2,2,1), k=c(3,3,2),
                     x=runif(3), dims=rep(3, 3))
  #' cumulate entries in duplicate positions
  Ms <- sparse3Darray(i=c(1,1,2), j=c(2,2,1), k=c(3,3,2),
                      x=runif(3), dims=rep(3, 3), strict=TRUE)

  #' print method
  print(M)
  
  #' conversion of other data
  A <- array(c(1,3,0,0,0,0,0,4,0,2,0,5,
               0,0,1,0,0,0,1,0,0,0,1,0),
             dim=c(3,4,2))
  A1 <- A[,,1]
  A2 <- A[,,2]
  Z <- A[integer(0), , ]
  
  #' array to sparse array
  AA <- as.sparse3Darray(A) # positive extent
  ZZ <- as.sparse3Darray(Z) # zero extent
  #' list of matrices to sparse array
  AA <- as.sparse3Darray(list(A1, A2))
  #' matrix to sparse array
  AA1 <- as.sparse3Darray(A1)
  #' vector to sparse array
  A11 <- A[,1,1]
  AA11 <- as.sparse3Darray(A11)
  #' NULL with warning
  as.sparse3Darray(list())

  #' 
  dim(AA) <- dim(AA) + 1

  I1 <- SparseIndices(A1)
  I11 <- SparseIndices(A11)
  
  if(require(Matrix)) {
    #' sparse matrices from Matrix package
    A1 <- as(A1, "sparseMatrix")
    A2 <- as(A2, "sparseMatrix")
    A11 <- as(A11, "sparseVector")
    #' convert a list of sparse matrices to sparse array
    AA <- as.sparse3Darray(list(A1, A2))
    #' sparse matrix to sparse array
    AA1 <- as.sparse3Darray(A1)
    #' sparse vector to sparse array
    AA11 <- as.sparse3Darray(A11)

    #' internals 
    E1  <- SparseEntries(A1)
    I1  <- SparseIndices(A1)
    I11 <- SparseIndices(A11)
    df <- data.frame(i=c(1,3,5), j=3:1, k=rep(2, 3), x=runif(3))
    aa <- EntriesToSparse(df, NULL)
    bb <- EntriesToSparse(df, 7)
    cc <- EntriesToSparse(df, c(7, 4))
    dd <- EntriesToSparse(df, c(7, 4, 3))
    #' duplicated entries
    dfdup <- df[c(1:3, 2), ]
    aa <- EntriesToSparse(dfdup, NULL)
    bb <- EntriesToSparse(dfdup, 7)
    cc <- EntriesToSparse(dfdup, c(7, 4))
    dd <- EntriesToSparse(dfdup, c(7, 4, 3))
  }

  BB <- evalSparse3Dentrywise(AA + AA/2)

  MM <- bind.sparse3Darray(M, M, along=1)
  MM <- bind.sparse3Darray(M, M, along=2)

  RelevantEmpty(42)
})

    
local({

  if(require(Matrix)) {

    M <- sparse3Darray(i=1:4, j=sample(1:4, replace=TRUE),
                       k=c(1,2,1,2), x=1:4, dims=c(5,5,2))

    M

    dimnames(M) <- list(letters[1:5], LETTERS[1:5], c("yes", "no"))
    M
    
    U <- aperm(M, c(1,3,2))
    U

    #' tests of [.sparse3Darray
    M[ 3:4, , ]
    M[ 3:4, 2:4, ]
    M[ 4:3, 4:2, 1:2]
    M[, 3, ]
    M[, 3, , drop=FALSE]
    M[c(FALSE,TRUE,FALSE,FALSE,TRUE), , ]
    M[, , c(FALSE,FALSE), drop=FALSE]
    M[1:2, 1, 2:3] # exceeds array bounds
    # matrix index
    M[cbind(3:5, 3:5, c(1,2,1))]
    M[cbind(3:5, 3:5, 2)]
    M[cbind(3:5,   2, 2)]
    M[cbind(c(2,2,4), c(3,3,2), 1)] # repeated indices
    M[cbind(1:4, 1, 2:3)] # exceeds array bounds

    MA <- as.array(M)
    UA <- as.array(U)

    Mfix <- sparse3Darray(i=1:3, j=c(3,1,2), k=4:2,
                          x=runif(3), dims=rep(4, 3))
    Mfix[cbind(1,3,4)] # single entry - occupied
    Mfix[cbind(1,2,4)] # single entry - unoccupied
    Mfix[cbind(1,c(2,3,2,3),4)] # sparse vector with repeated entries
    

    ## tests of "[<-.sparse3Darray"
    Mflip <- Mzero <- MandM <- Mnew <- Mext <- M
    Mflip[ , , 2:1] <- M
    stopifnot(Mflip[3,1,1] == M[3,1,2])
    Mzero[1:3,1:3,] <- 0
    stopifnot(all(Mzero[1,1,] == 0))
    M2a <- M[,,2,drop=FALSE]
    M2d <- M[,,2,drop=TRUE]
    MandM[,,1] <- M2a
    MandM[,,1] <- M2d
    ## slices of different dimensions
    M[ , 3, 1] <- 1:5
    M[2,  , 2] <- 1:5
    M[ 1, 3:5, 2] <- 4:6
    M[ 2, 5:3, 2] <- 4:6
    V3 <- sparseVector(x=1, i=2, length=3)
    M[ 1, 3:5, 2] <- V3
    M[ 2, 5:3, 2] <- V3
    M[,,2] <- M2a
    M[,,2] <- (M2a + 1)
    V5 <- sparseVector(x=1:2, i=2:3, length=5)
    M[,2,2] <- V5
    M[,,2] <- V5
    Mext[1,2,3] <- 4 # exceeds array bounds
    ## integer matrix index
    Mnew[cbind(3:5, 3:5, c(1,2,1))] <- 1:3
    Mnew[cbind(3:5, 3:5, 2)] <- 1:3
    Mnew[cbind(3:5,   2, 2)] <- 1:3
    Mnew[cbind(3:5, 3:5, c(1,2,1))] <- V3
    Mnew[cbind(3:5, 3:5, 2)] <- V3
    Mnew[cbind(3:5,   2, 2)] <- V3
    ## tests of arithmetic (Math, Ops, Summary)
    negM <- -M
    oneM <- 1 * M
    oneM <- M * 1
    twoM <- M + M
    range(M)

    cosM <- cos(M)  # non-sparse
    sinM <- sin(M)  # sparse

    Mpos <- (M > 0) # sparse
    Mzero <- !Mpos # non-sparse

    stopifnot(all((M+M) == 2*M))     # non-sparse
    stopifnot(!any((M+M) != 2*M))    # sparse

    ztimesM <- (1:5) * M  # sparse
    zplusM <- (1:5) + M  # non-sparse

    ## reconcile dimensions
    Msub <- M[,,1,drop=FALSE]
    Mdif <- M - Msub
    Mduf <- Msub - M
    
    ## tensor operator
    o <- tensorSparse(c(1,-1), M, 1, 3)
    o <- tensorSparse(M, M, 1:2, 1:2)
    o <- tensorSparse(M, M, 1:2, 2:1)
    o <- tensorSparse(as.array(M), as.array(M), 1:2, 2:1)
    V <- sparseVector(i=c(1,3,6),x=1:3, length=7)
    o <- tensorSparse(V,V)
    o <- tensorSparse(V,V,1,1)
    o <- tensorSparse(M,V[1:5],1,1)
    A <- sparseMatrix(i=integer(0), j=integer(0), x=numeric(0), dims=c(7, 15))
    A[1:4, 2:5] <- 3
    o <- tensorSparse(A, A, 1, 1)
    o <- tensorSparse(t(A), A, 2, 1)
    o <- tensorSparse(V, A, 1, 1)
    o <- tensorSparse(t(A), V, 2, 1)
    o <- tensorSparse(as.vector(V), A, 1, 1)
    o <- tensorSparse(t(A), as.vector(V), 2, 1)

    v <- 0:3
    o <- tensor1x1(v, Mfix)
    o <- tensor1x1(v, as.array(Mfix))
    o <- tensor1x1(as(v, "sparseVector"), Mfix)
    
    ## test of anyNA method
    anyNA(M)

    ## previously caused an error 
    a <- list(i = c(1L, 1L, 1L, 1L, 1L, 1L, 1L,
                    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                    2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
                    2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
                    2L, 2L, 2L, 2L),
              j = c(17L, 4L, 34L, 39L, 38L, 25L, 14L, 
                    40L, 1L, 19L, 36L, 9L, 16L, 23L,
                    15L, 17L, 4L, 34L, 39L, 38L, 
                    25L, 14L, 40L, 1L, 19L, 36L, 9L,
                    16L, 23L, 15L, 13L, 31L, 8L, 
                    5L, 42L),
              k = c(14L, 8L, 38L, 30L, 17L, 5L, 9L,
                    6L, 31L, 39L, 26L, 27L, 41L, 1L,
                    28L, 14L, 8L, 38L, 30L, 17L, 5L, 9L, 6L, 31L, 
                    39L, 26L, 27L, 41L, 1L, 28L, 36L, 15L, 19L, 21L, 42L))
    A <- with(a, sparse3Darray(i=i, j=j, k=k, x=1, dims=c(2, 42, 42)))
    stopifnot(all(sumsymouterSparse(A) == sumsymouter(as.array(A))))
    
    # no entries indexed
    A[integer(0), integer(0), integer(0)] <- 99
    A[matrix(, 0, 3)] <- 99

    if(FULLTEST) { # re-check with randomised data 
      ## .......... a possible application in spatstat
      ## n <- npoints(cells)
      ## cl10 <- as.data.frame(closepairs(cells, 0.1))
      ## cl12 <- as.data.frame(closepairs(cells, 0.12))
      ## ...........
      n <- 42
      ii <- sample(1:n, 20)
      jj <- sample(1:n, 20)
      cl12 <- data.frame(i=ii, j=jj)
      cl10 <- data.frame(i=ii[1:15], j=jj[1:15])
      ## ...........
      cl10$k <- 1
      cl12$k <- 2
      cl <- rbind(cl10, cl12)
      Z <- with(cl, sparse3Darray(i=i, j=j, k=k, x=1, dims=c(n,n,2)))
      dimnames(Z) <- list(NULL, NULL, c("r=0.1", "r=0.12"))
      Z <- aperm(Z, c(3,1,2))
      stopifnot(all(sumsymouterSparse(Z) == sumsymouter(as.array(Z))))
    }
    
    ## complex valued arrays
    Mcplx <- sparse3Darray(i=1:3, j=c(3,1,2), k=4:2,
                           x=runif(3)+runif(3)*1i, dims=rep(4, 3))
    print(Mcplx)
    

    #' -----------  sparsecommon.R -----------------------
    B <- sparseMatrix(i=1:3, j=3:1, x= 10 * (1:3), dims=c(4,4))
    #' (and using sparse 3D array M and sparse vector V from above)
    V2 <- sparseVector(i=c(2,3,6),x=4:6, length=7)  # different pattern
    check.anySparseVector(V2, 10, fatal=FALSE)

    Bmap <- mapSparseEntries(B, 1, 4:1)
    Mmap1 <- mapSparseEntries(M, 1, 5:1, across=3)
    Mmap2 <- mapSparseEntries(M, 3, 2:1, conform=FALSE)
    Mmap3 <- mapSparseEntries(M, 1, matrix(1:10, 5, 2), across=3)
    
    Vmap <- mapSparseEntries(V, 1, V2)
    Vmap <- mapSparseEntries(V, 1, 8)
    Vthrice  <- expandSparse(V, 3)
    VthriceT <- expandSparse(V, 3, 1)
    VF <- as.vector(V) # non-sparse
    VFmap <- mapSparseEntries(VF, 1, V2)
    VFmap <- mapSparseEntries(VF, 1, 8)
    VFthrice  <- expandSparse(VF, 3)
    VFthriceT <- expandSparse(VF, 3, 1)
    VFthriceX <- expandSparse(VF, 3, 2)
    
    VV <- sparseVectorCumul(rep(1:3,2), rep(c(3,1,2), 2), 5)

    Vsum <- applySparseEntries(V, sum)
    Bdouble <- applySparseEntries(B, function(x) { 2 * x })
    Mminus <- applySparseEntries(M, function(x) -x)

    VX <- expandSparse(B, 3, 1)
    VX <- expandSparse(B, 3, 2)
    VX <- expandSparse(B, 3, 3)
    
    # empty sparse matrices/arrays
    Bempty <- B
    Bempty[] <- 0
    mapSparseEntries(Bempty, 1, 42)
    Mempty <- M
    Mempty[] <- 0
    Mmap1 <- mapSparseEntries(Mempty, 1, 5:1, across=3)
    Mmap2 <- mapSparseEntries(Mempty, 3, 2:1, conform=FALSE)
    Mmap3 <- mapSparseEntries(Mempty, 1, matrix(1:10, 5, 2), across=3)

    #'  -------------- sparselinalg.R -------------------------
    U <- aperm(M,c(3,1,2))  # 2 x 5 x 5
    UU <- sumsymouterSparse(U, dbg=TRUE)
    w <- matrix(0, 5, 5)
    w[cbind(1:3,2:4)] <- 0.5
    w <- as(w, "sparseMatrix")
    UU <- sumsymouterSparse(U, w, dbg=TRUE)
    Uempty <- sparse3Darray(dims=c(2,5,5))
    UU <- sumsymouterSparse(Uempty, w, dbg=TRUE)
    #' complex
    Ucom <- U + U * 1i
    UU <- sumsymouterSparse(Ucom)
    UU <- sumsymouterSparse(Ucom, w)
    #' 
  }

  ## 1 x 1 x 1 arrays
  M1 <- sparse3Darray(i=1, j=1, k=1, x=42, dims=rep(1,3))
  M0 <- sparse3Darray(                     dims=rep(1,3))
  i1 <- matrix(1, 1, 3)
  a1 <- M1[i1]
  a0 <- M0[i1]
  A <- array(runif(75) * (runif(75) < 0.7), dim=c(3,5,5))
  M <- as.sparse3Darray(A)
  M[rep(1,3), c(1,1,2), rep(2, 3)]
})

}
