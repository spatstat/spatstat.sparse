# spatstat.sparse

## Sparse 3-dimensional arrays and linear algebra utilities

[![CRAN status](http://www.r-pkg.org/badges/version/spatstat.sparse)](http://CRAN.R-project.org/package=spatstat.sparse)
[![GitHub R package version](https://img.shields.io/github/r-package/v/spatstat/spatstat.sparse)](https://github.com/spatstat/spatstat.sparse)
![R-CMD-check](https://github.com/spatstat/spatstat.sparse/workflows/R-CMD-check/badge.svg)
[![Code Coverage Score](https://codecov.io/github/spatstat/spatstat.sparse/coverage.svg?branch=master)](https://codecov.io/github/spatstat/spatstat.sparse?branch=master)

This repository contains the current _development version_ of the
`spatstat.sparse` package.

For the most recent _official release_ of `spatstat.sparse`,
see the [CRAN page](https://CRAN.R-project.org/package=spatstat.sparse). 

### Functionality provided

The `spatstat.sparse` package 

  - defines a class of sparse three-dimensional arrays
    and supports standard operations on them.

  - provides utility functions for matrix computations
    that are common in statistics,
    such as quadratic forms.

It supports

   - creation of sparse arrays from raw data (numeric, integer, logical, or complex values)
   
   - conversion to/from other data types
   
   - array indexing, extraction of entries, assignment of new values
   
   - arithmetic and logical operations
   
   - tensor operations (generalising matrix multiplication)
   
   - permutation of array dimensions
   
   - binding of several arrays into a single array
   
   - printing of sparse arrays.

   - linear algebra operations including calculation of quadratic forms,
     bilinear forms, fractional powers of a matrix.

