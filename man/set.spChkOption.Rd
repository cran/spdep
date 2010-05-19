\name{set.spChkOption}
\alias{set.spChkOption}
\alias{get.spChkOption}
\alias{chkIDs}
\alias{spNamedVec}
\alias{set.VerboseOption}
\alias{get.VerboseOption}
\alias{set.ClusterOption}
\alias{get.ClusterOption}
\alias{set.ZeroPolicyOption}
\alias{get.ZeroPolicyOption}
\alias{set.rlecuyerSeedOption}
\alias{get.rlecuyerSeedOption}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Control checking of spatial object IDs}
\description{
  Provides support for checking the mutual integrity of spatial neighbour weights and spatial data; similar mechanisms are used for passing global verbose and zero.policy options, and for providing access to a running cluster for embarrassingly parallel tasks. 
}
\usage{
set.spChkOption(check)
get.spChkOption()
chkIDs(x, listw)
spNamedVec(var, data)
set.VerboseOption(check)
get.VerboseOption()
set.ZeroPolicyOption(check)
get.ZeroPolicyOption()
set.ClusterOption(cl)
get.ClusterOption()
set.rlecuyerSeedOption(seed)
get.rlecuyerSeedOption()
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{check}{a logical value, TRUE or FALSE}
  \item{cl}{a cluster object created by \code{makeCluster} in \pkg{snow}}
  \item{seed}{an integer vector of length 6 to seed the Lecuyer streams for cluster RNG}
  \item{x}{a vector the same length, or a two-dimensional array, or data frame with the same number of rows as the neighbours list in listw}
  \item{listw}{a \code{listw} object or \code{nb} object inheriting from "nb"}
  \item{var}{a character string or integer value for the column to be selected}
  \item{data}{a two-dimensional array or data frame containing var}
}
\details{
  Analysis functions will have an spChk argument by default set to NULL, and will call \code{get.spChkOption()} to get the global spatial option for whether to check or not --- this is initialised to FALSE, and consequently should not break anything. It can be changed to TRUE using \code{set.spChkOption(TRUE)}, or the spChk argument can be assigned in analysis functions. \code{spNamedVec()} is provided to ensure that rownames are passed on to single columns taken from two-dimensional arrays and data frames.
}
\value{
  \code{set.spChkOption()} returns the old logical value, \code{get.spChkOption()} returns the current logical value, and \code{chkIDs()} returns a logical value for the test lack of difference. \code{spNamedVec()} returns the selected column with the names set to the row names of the object from which it has been extracted.
}

\author{Roger Bivand \email{Roger.Bivand@nhh.no}}
\note{
The motivation for this mechanism is provided by the observation that spatial objects on a map and their attribute data values need to be linked uniquely, to avoid spurious results. The reordering between the legacy Columbus data set used the earlier publications and that available for download from the Spacestat website is just one example of a common problem.
}

 
\examples{
data(oldcol)
rownames(COL.OLD)
data(columbus)
rownames(columbus)
get.spChkOption()
oldChk <- set.spChkOption(TRUE)
get.spChkOption()
chkIDs(COL.OLD, nb2listw(COL.nb))
chkIDs(columbus, nb2listw(col.gal.nb))
chkIDs(columbus, nb2listw(COL.nb))
tmp <- try(moran.test(spNamedVec("CRIME", COL.OLD), nb2listw(COL.nb)))
print(tmp)
tmp <- try(moran.test(spNamedVec("CRIME", columbus), nb2listw(col.gal.nb)))
print(tmp)
tmp <- try(moran.test(spNamedVec("CRIME", columbus), nb2listw(COL.nb)))
print(tmp)
set.spChkOption(FALSE)
get.spChkOption()
moran.test(spNamedVec("CRIME", columbus), nb2listw(COL.nb))
tmp <- try(moran.test(spNamedVec("CRIME", columbus), nb2listw(COL.nb),
 spChk=TRUE))
print(tmp)
set.spChkOption(oldChk)
get.spChkOption()
}
\keyword{ spatial }
