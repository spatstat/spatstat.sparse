#'
#'  sparseMarkov.R
#'
#'  discrete-time finite-state Markov chain simulation
#'  using a sparse matrix representation of the transition matrix
#'
#'  $Revision: 1.5 $ $Date: 2026/04/11 04:25:09 $
#'
#'  Copyright (c) Adrian Baddeley 2026
#'  GNU Public Licence (>= 2.0)


runSparseMarkovChain <- function(P, x0, nsteps, ..., 
                                 result=c("history", "last"),
                                 check=TRUE,
                                 method=c("C", "interpreted")) {
  P <- as(P, "RsparseMatrix")
  if(!inherits(P, "RsparseMatrix"))
    stop("Unable to convert P to a sparse matrix in row-major form",
         call.=FALSE)
  method <- match.arg(method)
  result <- match.arg(result)
  savehistory <- (result == "history")
  x0 <- as.integer(x0)
  nx <- length(x0)
  if(check) {
    ra <- range(P)
    if(!all(is.finite(ra)))
      stop("P contains infinite, NA or NaN entries", call.=FALSE)
    if(ra[1L] < 0) stop("P contains negative entries", call.=FALSE)
    if(ra[2L] > 1) stop("P contains entries greater than 1", call.=FALSE)
    rs <- range(rowSums(P))
    if(max(abs(rs-1)) > sqrt(.Machine$double.eps))
      stop("Some of the rows of P do not sum to 1", call.=FALSE)
    rx <- range(x0)
    if(rx[1L] < 1 || rx[2L] > nrow(P))
      stop("Some indices in x0 are out of range", call.=FALSE)
  }
  ## >>>>>>>>>>>>  run chain  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  switch(method,
         interpreted = {
           x <- x0
           if(savehistory) history <- x0
           for(istep in 1:nsteps) {
             xnew <- integer(0)
             ## make one transition for each particle
             for(i in x) {
               p.i <- P[i, ]
               jj <- which(p.i > 0)
               p.ij <- p.i[jj]
               xnew <- c(xnew, sample(jj, 1L, prob=p.ij))
             }
             ## current state of each particle
             x <- xnew
             ## append to history (column = step, row = particle)
             if(savehistory) history <- cbind(history, x)
           }
         },
         C = {
           if(savehistory) {
             z <- .C(SP_rMCspMH,
                     nrows = as.integer(dim(P)[1L]),
                     nsparse = as.integer(length(P@x)),
                     probval = as.double(P@x),
                     colindex = as.integer(P@j),
                     rowstart = as.integer(P@p),
                     npoints = as.integer(nx),
                     startpos = as.integer(x0 - 1L),  # zero-based indices
                     nsteps = as.integer(nsteps),
                     history = as.integer(integer(nx * (nsteps + 1L))),
                     PACKAGE = "spatstat.sparse")
             history <- matrix(z$history + 1L, nx, nsteps+1L, byrow=TRUE)
           } else {
             z <- .C(SP_rMCspMF,
                     nrows = as.integer(dim(P)[1L]),
                     nsparse = as.integer(length(P@x)),
                     probval = as.double(P@x),
                     colindex = as.integer(P@j),
                     rowstart = as.integer(P@p),
                     npoints = as.integer(nx),
                     startpos = as.integer(x0 - 1L), # zero-based indices
                     nsteps = as.integer(nsteps),
                     endpos = as.integer(integer(nx)),
                     PACKAGE = "spatstat.sparse")
             x <- z$endpos + 1L
           }
         })
  if(savehistory) {
    dimnames(history) <- NULL
    if(nx == 1) history <- history[1L, ]
    return(history)
  }
  return(x)
}


  
