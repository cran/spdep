% Copyright 2003 by Roger S. Bivand
\name{EBImoran.mc}
\alias{EBImoran.mc}
\alias{EBImoran}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Permutation test for empirical Bayes index}
\description{
  An empirical Bayes index modification of Moran's I for testing for spatial autocorrelation in a rate, typically the number of observed cases in a population at risk. The index value is tested by using nsim random permutations of the index for the given spatial weighting scheme, to establish the rank of the observed statistic in relation to the nsim simulated values.
}
\usage{
EBImoran.mc(n, x, listw, nsim, zero.policy = FALSE, 
 alternative = "greater", spChk=NULL)
}
EBImoran(z, listw, nn, S0, zero.policy=FALSE)
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{n}{a numeric vector of counts of cases the same length as the neighbours list in listw}
  \item{x}{a numeric vector of populations at risk the same length as the neighbours list in listw}
  \item{listw}{a \code{listw} object created for example by \code{nb2listw}}
  \item{nsim}{number of permutations}
  \item{zero.policy}{if TRUE assign zero to the lagged value of zones without neighbours, if FALSE assign NA}
  \item{alternative}{a character string specifying the alternative hypothesis, must be one of greater (default), less or two.sided}
  \item{spChk}{should the data vector names be checked against the spatial objects for identity integrity, TRUE, or FALSE, default NULL to use \code{get.spChkOption()}}}
  \item{z}{a vector of indices based on n and x}
  \item{nn}{the number of observations}
  \item{S0}{global sum of weights}
}
\details{
The statistic used is (m is the number of observations):
\deqn{EBI = \frac{m}{\sum_{i=1}^{m}\sum_{j=1}^{m}w_{ij}}
\frac{\sum_{i=1}^{m}\sum_{j=1}^{m}w_{ij}z_i z_j}{\sum_{i=1}^{m}(z_i - \bar{z})^2}
}{EBI = (n sum_i sum_j w_ij z_i z_ji) / (S0 sum_i (z_i - zbar)^2)}
where:
\deqn{z_i = \frac{p_i - b}{\sqrt{v_i}}}{z_i = (p_i - b) / sqrt(v_i)}
and:
\deqn{p_i = n_i / x_i}{p_i = n_i / x_i}
\deqn{v_i = a + (b / x_i)}{v_i = a + (b / x_i)}
\deqn{b = \sum_{i=1}^{m} n_i / \sum_{i=1}^{m} x_i }{sum_i n_i / sum_i x_i}
\deqn{a = s^2 - b / (\sum_{i=1}^{m} x_i / m)}{a = s^2 - b / (sum_i x_i / m)}
\deqn{s^2 = \sum_{i=1}^{m} x_i (p_i - b)^2 / \sum_{i=1}^{m} x_i }{s^2 = sum_i x_i (p_i - b)^2 / sum_i x_i}
  
}
\value{
A list with class \code{htest} and \code{mc.sim} containing the 
following components:
  \item{statistic}{the value of the observed Moran's I.}
  \item{parameter}{the rank of the observed Moran's I.}
  \item{p.value}{the pseudo p-value of the test.}
  \item{alternative}{a character string describing the alternative hypothesis.}
  \item{method}{a character string giving the method used.}
  \item{data.name}{a character string giving the name(s) of the data, and the number of simulations.}
  \item{res}{nsim simulated values of statistic, final value is observed statistic}
  \item{z}{a numerical vector of Empirical Bayes indices as z above}

}
\references{Assunção RM, Reis EA 1999 A new proposal to adjust Moran's I for population density. Statistics in Medicine 18, pp. 2147--2162}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{moran}}, \code{\link{moran.mc}}, \code{\link{BH}}}

\examples{
data(nc.sids)
EBImoran.mc(spNamedVec("SID74", nc.sids), spNamedVec("BIR74", nc.sids),
 nb2listw(sidsorig.nb, style="B", zero.policy=TRUE), nsim=999,
 zero.policy=TRUE)
sids.p <- nc.sids$SID74 / nc.sids$BIR74
names(sids.p) <- rownames(nc.sids)
moran.mc(sids.p, nb2listw(sidsorig.nb, style="B", zero.policy=TRUE),
 nsim=999, zero.policy=TRUE)
}
\keyword{spatial}% at least one, from doc/KEYWORDS

