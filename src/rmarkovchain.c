/*
  rmarkovchain.c

  Simulate a finite-state discrete-time Markov chain
  using sparse representation of transition matrix 

  rMCspMF  - return final state only
  rMCspMH  - return history of chain

  These C functions are called from 'sparseMarkov.R'

  The transition matrix M is required to be in row-major sparse form
  (virtual class 'RsparseMatrix' which embraces 'dgRMatrix' and 'dtRMatrix')
  Matrices can be converted to this form using 'as(M, "RsparseMatrix")'

  Data passed to the C code from the sparse matrix M:
  
  nrows     Number of rows (and columns) of matrix            dim(M)[1]
  nsparse   Number of nonzero entries in matrix               length(M@x)
  probval   Vector of nonzero entries                         M@x
  colindex  Vector of column indices for nonzero entries      M@j
  rowstart  Start position in M@x, M@j of each row of matrix  M@p

  The indices in 'colindex' and 'rowstart' are zero-based.

  $Revision: 1.8 $ $Date: 2026/04/11 09:54:48 $

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
  register int steps, points, entries, currentpos, j;
  int thisrowstart, nextrowstart, thisrowlength;
  register double u;
  int Nrows, Nsparse, Npoints, Nsteps;
  register int *sp, *ep, *cp;
  register double *pp;

  Nrows   = *nrows;
  Nsparse = *nsparse;
  Npoints = *npoints;
  Nsteps  = *nsteps;

  GetRNGstate();

  /* initialise pointers */
  sp = startpos;
  ep = endpos;
  
  for(points = Npoints; points > 0; --points, ++sp, ++ep) {
    /* initialise position */
    currentpos = *sp;
    /* run chain */
    for(steps = Nsteps; steps > 0; --steps) {
      /* nonzero transition probabilities from current position */
      thisrowstart = rowstart[currentpos];
      nextrowstart = ((currentpos+1) < Nrows) ? rowstart[currentpos+1] : Nsparse;
      thisrowlength = nextrowstart - thisrowstart;
      /* random number */
      u = unif_rand();
      /* pointers into probval[] and colindex[] for this row */
      pp = probval  + thisrowstart;
      cp = colindex + thisrowstart;
      /* loop over nonzero entries in current row */
      j = -1;
      for(entries = thisrowlength; entries > 0; --entries, pp++, cp++) {
	u -= *pp;
	if(u <= 0.0) {
	  /* jump */
	  j = *cp;
	  break;
	}
      }
      if(j < 0) {
	/* 
	   Random number exceeded row sum of transition matrix.
	   theoretically impossible -- can occur due to numerical error 
	*/
	j = colindex[nextrowstart - 1];
      }
      currentpos = j;
    }
    /* save state */
    *ep = currentpos;
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
  register int steps, points, entries, currentpos, j;
  int thisrowstart, nextrowstart, thisrowlength;
  register double u;
  int Nrows, Nsparse, Npoints, Nsteps, Nhistory;
  register int *sp, *cp, *histp;
  register double *pp;
  
  Nrows   = *nrows;
  Nsparse = *nsparse;
  Npoints = *npoints;
  Nsteps  = *nsteps;

  Nhistory = Nsteps + 1;
  
  GetRNGstate();

  /* pointer to next entry in 'history' */
  histp = history;
  /* pointer to initial state */
  sp = startpos;
  
  for(points = Npoints; points > 0; --points, ++sp) {
    /* initialise position */
    currentpos = *sp;
    /* save initial state */
    *histp = currentpos;
    ++histp;
    /* run chain */
    for(steps = Nsteps; steps > 0; --steps, ++histp) {
      /* nonzero transition probabilities from current position */
      thisrowstart = rowstart[currentpos];
      nextrowstart = ((currentpos+1) < Nrows) ? rowstart[currentpos+1] : Nsparse;
      thisrowlength = nextrowstart - thisrowstart;
      /* random number */
      u = unif_rand();
      /* pointers into probval[] and colindex[] for this row */
      pp = probval  + thisrowstart;
      cp = colindex + thisrowstart;
      /* loop over nonzero entries in current row */
      j = -1;
      for(entries = thisrowlength; entries > 0; --entries, pp++, cp++) {
	u -= *pp;
	if(u <= 0.0) {
	  /* jump */
	  j = *cp;
	  break;
	}
      }
      if(j < 0) {
	/* 
	   Random number exceeded row sum of transition matrix.
	   theoretically impossible -- can occur due to numerical error 
	*/
	j = colindex[nextrowstart - 1];
      }
      currentpos = j;
      /* save state */
      *histp = currentpos;
    }
  }
  
  PutRNGstate();
}

	      
