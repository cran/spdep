% Copyright 2002 by Roger S. Bivand
\name{sp.correlogram}
\alias{sp.correlogram}
\alias{plot.spcor}
\alias{print.spcor}
\title{Spatial correlogram}
\description{
  Spatial correlograms for Moran's I and the autocorrelation coefficient, with
print and plot helper functions.
}
\usage{
sp.correlogram(neighbours, var, order = 1, method = "corr", style = "W", 
  randomisation = TRUE, zero.policy = FALSE, spChk=NULL)
plot.spcor(x, main, ylab, ylim, ...)
print.spcor(x, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{neighbours}{an object of class \code{nb}}
  \item{var}{a numeric vector}
  \item{order}{maximum lag order}
  \item{method}{"corr" for correlation, "I" for Moran's I}
  \item{style}{\code{style} can take values W, B, C, and S}
  \item{randomisation}{variance of I calculated under the assumption of randomisation, if FALSE normality}
  \item{zero.policy}{If FALSE stop with error for any empty neighbour sets, 
if TRUE permit the weights list to be formed with zero-length weights vectors}
  \item{spChk}{should the data vector names be checked against the spatial objects for identity integrity, TRUE, or FALSE, default NULL to use \code{get.spChkOption()}}
  \item{x}{an object from \code{sp.correlogram()} of class \code{spcor}}
  \item{main}{an overall title for the plot}
  \item{ylab}{a title for the y axis}
  \item{ylim}{the y limits of the plot}
  \item{\dots}{further arguments passed through}
}


\value{
  returns a list of class \code{spcor}:
  \item{res}{for "corr" a vector of values; for "I", a matrix of estimates of "I", expectations, and variances}
  \item{method}{"I" or "corr"}
  \item{cardnos}{list of tables of neighbour cardinalities for the lag orders used}
  \item{var}{variable name}
}
\references{Cliff, A. D., Ord, J. K. 1981 \emph{Spatial processes}, Pion, pp. 
118--122, Martin, R. L., Oeppen, J. E. 1975 The identification of regional
forecasting models using space-time correlation functions, \emph{Transactions
of the Institute of British Geographers}, 66, 95--118.}
\author{Roger Bivand, \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{nblag}}, \code{\link{moran}}}

\examples{
data(nc.sids)
ft.SID74 <- sqrt(1000)*(sqrt(nc.sids$SID74/nc.sids$BIR74) +
  sqrt((nc.sids$SID74+1)/nc.sids$BIR74))
tr.SIDS74 <- ft.SID74*sqrt(nc.sids$BIR74)
names(tr.SIDS74) <- rownames(nc.sids)
print(sp.correlogram(sidsorig.nb, tr.SIDS74, order=8, method="corr",
 zero.policy=TRUE))
print(sp.correlogram(sidsorig.nb, tr.SIDS74, order=8, method="I",
 zero.policy=TRUE))
plot(sp.correlogram(sidsorig.nb, tr.SIDS74, order=8, method="I",
 zero.policy=TRUE))
plot(sp.correlogram(sidsorig.nb, tr.SIDS74, order=8, method="corr",
 zero.policy=TRUE))
drop.no.neighs <- !(1:length(sidsorig.nb) \%in\% which(card(sidsorig.nb) == 0))
sub.sidsorig.nb <- subset(sidsorig.nb, drop.no.neighs)
plot(sp.correlogram(sub.sidsorig.nb, subset(tr.SIDS74,  drop.no.neighs),
 order=8, method="corr"))
}
\keyword{spatial}
