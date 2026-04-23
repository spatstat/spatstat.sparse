#' Header for spatstat.sparse/tests/*R
#'

require(spatstat.sparse)
ALWAYS <- FULLTEST <- TRUE
#'    tests/sparsemarkov.R
#'  Tests of code for Markov chain with sparse transition matrix
#'  $Revision: 1.3 $ $Date: 2026/04/23 05:20:03 $

if(!exists("ALWAYS")) ALWAYS <- TRUE
if(!exists("FULLTEST")) FULLTEST <- ALWAYS

if(FULLTEST) { 
  local({
    testit <- function(P, Pname, xstart=50, ns=10, np=3,
                       walkies=FALSE, absorb=NULL) {
      cat(paste(">>>>>> ", Pname, " <<<<<<<\n"))
      cat("\tsingle particle, final state..\n")
      X1 <- runSparseMarkovChain(P, x0=xstart, nsteps=ns, result="l")
      cat("\tsingle particle, history..\n")
      X2 <- runSparseMarkovChain(P, x0=xstart, nsteps=ns, result="h")
      cat("\tseveral particles, final states..\n")
      X3 <- runSparseMarkovChain(P, x0=rep(xstart, np), nsteps=ns,
                                 result="l")
      cat("\tseveral particles, histories..\n")
      X4 <- runSparseMarkovChain(P, x0=rep(xstart, np), nsteps=ns,
                                 result="h")
      if(length(X1) != 1)
        stop(paste(Pname, "end state was not a single particle"))
      if(length(X2) != (ns + 1))
        stop(paste(Pname, "history has wrong length"))
      if(length(X3) != np)
        stop(paste(Pname, "final number of particles was not preserved"))
      if(!is.matrix(X4))
        stop(paste(Pname, "a matrix was expected for the history"))
      if(nrow(X4) != np)
        stop(paste(Pname, "number of particles was not preserved"))
      if(ncol(X4) != (ns + 1))
        stop(paste(Pname, "histories have wrong length"))

      if(walkies) {
        ## jumps are all +-1
        if(abs(X1 - xstart) > ns)
          stop(paste(Pname, "wandered impossibly far"))
        if(!all(abs(diff(X2)) == 1))
          stop(paste(Pname, "jumps were not all +- 1"))
        if(max(abs(X3 - xstart)) > ns)
          stop(paste(Pname, "some particles wandered impossibly far"))
        if(!all(abs(apply(X4, 1, diff)) == 1))
          stop(paste(Pname, "jumps of particles were not all +- 1"))
      }

      if(!is.null(absorb)) {
        if(X1 != absorb)
          stop("Absorbing chain didn't absorb (final)")
        if(X2[ns+1] != absorb)
          stop("Absorbing chain didn't absorb (history)")
        if(any(X3 != absorb))
          stop("Absorbing chains didn't absorb (final)")
        if(any(X4[,ns+1] != absorb))
          stop("Absorbing chains didn't absorb (history)")
      }
      
      cat("OK\n")
      return(invisible(list(X1=X1, X2=X2, X3=X3, X4=X4)))
    }
      
    #' Simple random walk
    Pwalk <- matrix(0, 100, 100)
    Pwalk[abs(row(Pwalk) - col(Pwalk)) == 1] <- 1
    Pwalk <- Pwalk/rowSums(Pwalk)
    testit(Pwalk, "Simple random walk", walkies=TRUE)

    #' Absorbing
    Pabsorb <- matrix(0, 10, 10)
    Pabsorb[row(Pabsorb) < col(Pabsorb)] <- 1
    Pabsorb[,10] <- 1
    Pabsorb <- Pabsorb/rowSums(Pabsorb)
    a <- testit(Pabsorb, "Absorbing chain", xstart=1, ns=10, absorb=10)
  }
)}
