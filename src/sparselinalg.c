#include <R.h>
#include <R_ext/Utils.h>

/*
  sparselinalg.c

  Counterpart of 'linalg.c' for sparse matrices/arrays

  $Revision: 1.8 $  $Date: 2022/10/22 10:09:51 $

  Copyright (C) Adrian Baddeley, Ege Rubak and Rolf Turner 2001-2018
  Licence: GNU Public Licence >= 2

 */

#undef DBG

#define FNAME CspaSumSymOut
#undef WEIGHTS
#include "spasumsymout.h"
#undef FNAME

#define FNAME CspaWtSumSymOut
#define WEIGHTS
#include "spasumsymout.h"
#undef FNAME

#define DBG

#define FNAME CDspaSumSymOut
#undef WEIGHTS
#include "spasumsymout.h"
#undef FNAME

#define FNAME CDspaWtSumSymOut
#define WEIGHTS
#include "spasumsymout.h"
#undef FNAME

