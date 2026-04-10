/*
  rmarkovchain.c

  Simulate a finite-state discrete-time Markov chain
  using sparse representation of transition matrix 

  rMCspMF  - return final state only
  rMCspMH  - return history of chain

  Sparse matrix is required to be in row-major sparse form
  (class 'dgRmatrix') which can be achieved using 'as(, "RsparseMatrix")'
  Information passed to the C code:
  
  nrows     Number of rows (and columns) of matrix         dim(M)[1]
  nsparse   Number of nonzero entries in matrix            length(M@x)
  probval   Vector of nonzero entries                      M@x
  colindex  Vector of column indices for nonzero entries   M@j
  rowstart  Start position in M@x, M@j of each row         M@p

  The indices in 'colindex' and 'rowstart' are zero-based.

  $Revision: 1.4 $ $Date: 2026/04/10 08:04:44 $

  Copyright (c) Adrian Baddeley 2026
  GNU Public Licence (>= 2.0)

*/

#include <R.h>
#include <R_ext/Utils.h>
#include <math.h>

/* rMCspMF
   Run 'M'ultiple independent realisations
   from different starting points, 
   and return only the 'F'inal states 
*/

void rMCspMF(
	      /* transition matrix P in sparse form (row-major) */
	      int *nrows,           /* dimensions of P */
	      int *nsparse,         /* number of nonzero entries in P */
	      double *probval,      /* nonzero entries */
	      int *colindex,        /* column index for each nonzero value */
	      int *rowstart,        /* index of start of each row */
	      /* initial state */
	      int *npoints,         
	      int *startpos,        /* initial positions (zero-based) */
	      /* number of steps */
	      int *nsteps,
	      /* output: final state for each point */
	      int *endpos) {
  register int istep, ipoint, currentpos, thisrowstart, nextrowstart, j, k;
  register double u;
  int Nrows, Nsparse, Npoints, Nsteps;

  Nrows   = *nrows;
  Nsparse = *nsparse;
  Npoints = *npoints;
  Nsteps  = *nsteps;

  GetRNGstate();
  
  for(ipoint = 0; ipoint < Npoints; ipoint++) {
    /* initialise position */
    currentpos = startpos[ipoint];
    /* run chain */
    for(istep = 0; istep < Nsteps; istep++) {
      /* transition probabilities from current position */
      thisrowstart = rowstart[currentpos];
      nextrowstart = ((currentpos+1) < Nrows) ? rowstart[currentpos+1] : Nsparse;
      /* random number */
      u = unif_rand();
      j = -1;
      for(k = thisrowstart; k < nextrowstart; k++) {
	u -= probval[k];
	if(u <= 0.0) {
	  /* jump */
	  j = colindex[k];
	  break;
	}
      }
      if(j < 0) j = colindex[nextrowstart - 1];
      currentpos = j;
    }
    /* save state */
    endpos[ipoint] = currentpos;
  }
  
  PutRNGstate();
}

	      
/* rMCspMH
   Run 'M'ultiple independent realisations
   from different starting points, 
   and return the 'H'istory of each point.
*/

void rMCspMH(
	      /* transition matrix P in sparse form (row-major) */
	      int *nrows,           /* dimensions of P */
	      int *nsparse,         /* number of nonzero entries in P */
	      double *probval,      /* nonzero entries */
	      int *colindex,        /* column index for each nonzero value */
	      int *rowstart,        /* index of start of each row */
	      /* initial state */
	      int *npoints,         
	      int *startpos,        /* initial positions (zero-based) */
	      /* number of steps */
	      int *nsteps,
	      /* output: history for each point ( (Nsteps + 1) * Npoints ) */
	      int *history) {
  register int istep, ipoint, currentpos, thisrowstart, nextrowstart, j, k;
  register double u;
  int Nrows, Nsparse, Npoints, Nsteps, Nhistory;
  int *histp;

  Nrows   = *nrows;
  Nsparse = *nsparse;
  Npoints = *npoints;
  Nsteps  = *nsteps;

  Nhistory = Nsteps + 1;
  
  GetRNGstate();

  /* pointer to next entry in 'history' */
  histp = history;
  
  for(ipoint = 0; ipoint < Npoints; ipoint++) {
    /* initialise position */
    currentpos = startpos[ipoint];
    /* save initial state */
    *histp = currentpos;
    ++histp;
    /* run chain for point number 'ipoint' */
    for(istep = 0; istep < Nsteps; istep++, histp++) {
      /* transition probabilities from current position */
      thisrowstart = rowstart[currentpos];
      nextrowstart = ((currentpos+1) < Nrows) ? rowstart[currentpos+1] : Nsparse;
      /* random number */
      u = unif_rand();
      j = -1;
      for(k = thisrowstart; k < nextrowstart; k++) {
	u -= probval[k];
	if(u <= 0.0) {
	  /* jump */
	  j = colindex[k];
	  break;
	}
      }
      if(j < 0) j = colindex[nextrowstart - 1];
      currentpos = j;
      /* save state */
      *histp = currentpos;
    }
  }
  
  PutRNGstate();
}

	      
