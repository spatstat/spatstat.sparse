\name{spatstat.sparse-internal} 
\title{Internal Functions of spatstat.sparse Package}
\alias{spatstat.sparse-internal} %DoNotExport
%% indices
\alias{representativeRows}
%%%%  undocumented linear algebra
\alias{checksolve}
\alias{check.mat.mul}
\alias{sumsymouter}
%% sparse 3D arrays
\alias{unionOfSparseIndices}
\alias{inside3Darray}
\alias{SparseEntries}
\alias{SparseIndices}
\alias{EntriesToSparse}
\alias{mapSparseEntries}
\alias{applySparseEntries}
\alias{sumsymouterSparse}
\alias{rbindCompatibleDataFrames}
\alias{check.anySparseVector}
\alias{evalSparse3Dentrywise}
\alias{expandSparse}
\alias{sparseVectorCumul}
\alias{tensor1x1}
\description{
  Internal utility functions of the \code{spatstat.sparse} package.
}
\usage{
%% indices
representativeRows(x)
%%% undocumented linear algebra
checksolve(M, action, descrip, target)
check.mat.mul(A, B, Acols, Brows, fatal)
sumsymouter(x, w, distinct)
%% sparse 3D arrays
unionOfSparseIndices(A,B)
inside3Darray(d, i, j, k)
SparseEntries(x)
SparseIndices(x)
EntriesToSparse(df, dims)
mapSparseEntries(x, margin, values, conform, across)
applySparseEntries(x, f, \dots)
sumsymouterSparse(x, w, distinct, dbg)
rbindCompatibleDataFrames(x)
check.anySparseVector(v, npoints, fatal, things,
                      naok, warn, vname, oneok)
evalSparse3Dentrywise(expr, envir)
expandSparse(x, n, across)
sparseVectorCumul(x, i, length)
tensor1x1(A, B)
}
\details{
  These internal \pkg{spatstat.sparse} functions are not usually called
  directly by the user. Their names and capabilities may change
  without warning from one version of \pkg{spatstat.utils} to the next.
}
\keyword{internal}
