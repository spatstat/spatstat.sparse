\name{spatstat.sparse-internal} 
\title{Internal Functions of spatstat.sparse Package}
%% indices
\alias{fullIndexSequence}
\alias{grokIndexVector}
\alias{replacementIndex}
\alias{representativeRows}
\alias{strictIndexSequence}
%%%%  undocumented linear algebra
\alias{checksolve}
\alias{check.mat.mul}
\alias{sumsymouter}
%% sparse 3D arrays
\alias{dim.sparse3Darray}
\alias{dim<-.sparse3Darray}
\alias{dimnames.sparse3Darray}
\alias{dimnames<-.sparse3Darray}
\alias{print.sparse3Darray}
\alias{as.array.sparse3Darray}
\alias{[.sparse3Darray}
\alias{[<-.sparse3Darray}
\alias{anyNA.sparse3Darray}
\alias{RelevantZero}
\alias{RelevantNA}
\alias{RelevantEmpty}
\alias{isRelevantZero}
\alias{unionOfSparseIndices}
\alias{inside3Darray}
\alias{SparseEntries}
\alias{SparseIndices}
\alias{EntriesToSparse}
\alias{mapSparseEntries}
\alias{applySparseEntries}
\alias{sumsymouterSparse}
\alias{marginSumsSparse}
\alias{rbindCompatibleDataFrames}
\alias{bind.sparse3Darray}
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
fullIndexSequence(g)
grokIndexVector(ind, len, nama)
replacementIndex(ii, stuff)
representativeRows(x)
strictIndexSequence(g)
%%% undocumented linear algebra
checksolve(M, action, descrip, target)
check.mat.mul(A, B, Acols, Brows, fatal)
sumsymouter(x, w, distinct)
%% sparse 3D arrays
\method{dim}{sparse3Darray}(x)
\method{dim}{sparse3Darray}(x) <- value
\method{dimnames}{sparse3Darray}(x)
\method{dimnames}{sparse3Darray}(x) <- value
\method{print}{sparse3Darray}(x, \dots)
\method{as.array}{sparse3Darray}(x, \dots)
\method{[}{sparse3Darray}(x, i, j, k, drop, \dots)
\method{[}{sparse3Darray}(x, i, j, k, \dots) <- value
\method{anyNA}{sparse3Darray}(x, recursive)
RelevantZero(x)
RelevantNA(x)
RelevantEmpty(x)
isRelevantZero(x)
unionOfSparseIndices(A,B)
inside3Darray(d, i, j, k)
SparseEntries(x)
SparseIndices(x)
EntriesToSparse(df, dims)
mapSparseEntries(x, margin, values, conform, across)
applySparseEntries(x, f, \dots)
sumsymouterSparse(x, w, distinct, dbg)
marginSumsSparse(X, MARGIN)
rbindCompatibleDataFrames(x)
bind.sparse3Darray(A, B, along)
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
