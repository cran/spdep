% Copyright 2003 by Roger S. Bivand
\name{EBlocal}
\alias{EBlocal}
\title{Local Empirical Bayes estimator}
\description{
  The function computes local empirical Bayes estimates for rates "shrunk" to a neighbourhood mean for neighbourhoods given by the \code{nb} neighbourhood list.
}
\usage{
EBlocal(ri, ni, nb, zero.policy = FALSE, spChk = NULL, geoda=FALSE)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{ri}{a numeric vector of counts of cases the same length as the neighbours list in nb}
  \item{ni}{a numeric vector of populations at risk the same length as the neighbours list in nb}
  \item{nb}{a \code{nb} object of neighbour relationships}
  \item{zero.policy}{if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA}
  \item{spChk}{should the data vector names be checked against the spatial objects for identity integrity, TRUE, or FALSE, default NULL to use \code{get.spChkOption()}}
 \item{geoda}{default=FALSE, following Marshall's algorithm as interpreted by Bailey and Gatrell, pp. 305-307, and exercise 8.2, pp. 328-330 for the definition of phi; TRUE for the definition of phi used in GeoDa (see discussion on OpenSpace mailing list June 2003: http://agec221.agecon.uiuc.edu/pipermail/openspace/2003-June/thread.html)}
}
\details{
Details of the implementation are to be found in Marshall, p. 286, and Bailey and Gatrell p. 307 and exercise 8.2, pp. 328--330. The example results do not fully correspond to the sources because of slightly differing neighbourhoods, but are generally close.
  
}
\value{
A data frame with two columns:
  \item{raw}{a numerical vector of raw (crude) rates}
  \item{est}{a numerical vector of local empirical Bayes estimates}

and a \code{parameters} attribute list with components:

  \item{a}{a numerical vector of local phi values}
  \item{m}{a numerical vector of local gamma values}
}
\references{Marshall R M (1991) Mapping disease and mortality rates using Empirical Bayes Estimators, Applied Statistics, 40, 283--294; Bailey T, Gatrell A (1995) Interactive Spatial Data Analysis, Harlow: Longman, pp. 303--306.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}, based on contributions by Marilia Carvalho}

\seealso{\code{\link{EBest}}, \code{\link{probmap}}}

\examples{
data(auckland)
res <- EBlocal(spNamedVec("Deaths.1977.85", auckland),
 9*spNamedVec("Under.5.1981", auckland), auckland.nb)
brks <- c(-Inf,2,2.5,3,3.5,Inf)
cols <- grey(6:2/7)
library(maptools)
plot(auckpolys, col=cols[findInterval(res$est*1000, brks)], forcefill=FALSE)
legend(c(70,90), c(70,95), fill=cols, legend=leglabs(brks), bty="n")
title(main="Local moment estimator of infant mortality per 1000 per year")
}
\keyword{spatial}
