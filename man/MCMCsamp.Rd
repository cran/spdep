\name{MCMCsamp}
\alias{MCMCsamp}
\alias{MCMCsamp.spautolm}
\alias{MCMCsamp.sarlm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{MCMC sample from fitted spatial regression}
\description{The \code{MCMCsamp} method uses \code{\link[MCMCpack]{MCMCmetrop1R}} from \pkg{MCMCpack} to make MCMX samples from fitted maximum likelihood spatial regression models}
\usage{
MCMCsamp(object, mcmc = 1L, verbose = NULL, ...)
\method{MCMCsamp}{spautolm}(object, mcmc = 1L, verbose = NULL, ...,
 burnin = 0L, thin = 1L, seed = NA, verbose_step = 100L, listw,
 control = list())
\method{MCMCsamp}{sarlm}(object, mcmc = 1L, verbose = NULL, ...,
    burnin=0L, thin=1L, seed=NA, verbose_step=100L, listw, listw2=NULL,
    control=list())}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A spatial regression model object fitted by maximum likelihood with \code{\link{spautolm}}}
  \item{mcmc}{The number of MCMC iterations after burnin}
  \item{verbose}{default NULL, use global option value; if TRUE, reports progress}
  \item{\dots}{Arguments passed through}
  \item{burnin}{The number of burn-in iterations for the sampler}
  \item{thin}{The thinning interval used in the simulation.  The number of MCMC iterations must be divisible by this value}
  \item{seed}{default NA; if NA, the Mersenne Twister generator is used with default seed 12345; if an integer is passed it is used to seed the Mersenne twister}
  \item{verbose_step}{default 100; if \sQuote{verbose} is TRUE, the iteration number, the theta vector, the function value, and the Metropolis acceptance rate are sent to the screen every \sQuote{verbose}th iteration}
  \item{listw, listw2}{\code{listw} objects created for example by \code{nb2listw}; should be the same object(s) used for fitting the model}
  \item{control}{list of extra control arguments - see \code{\link{spautolm}}}
}
\value{An object of class \dQuote{mcmc} suited to \pkg{coda}}
\references{Martin AD, Quinn KM, Park JH (2011) MCMCpack: Markov Chain Monte Carlo in R, Journal of Statistical Software 42:9, 1-21, \url{http://www.jstatsoft.org/v42/i09}}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link[MCMCpack]{MCMCmetrop1R}}, \code{\link{spautolm}}}
\examples{
\dontrun{
example(NY_data)
esar1f <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY, family="SAR", method="eigen")
summary(esar1f)
res <- MCMCsamp(esar1f, mcmc=5000, burnin=500, listw=listw_NY)
summary(res)
ecar1f <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY, family="CAR", method="eigen")
summary(ecar1f)
res <- MCMCsamp(ecar1f, mcmc=5000, burnin=500, listw=listw_NY)
summary(res)
esar1fw <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY, weights=POP8, family="SAR", method="eigen")
summary(esar1fw)
res <- MCMCsamp(esar1fw, mcmc=5000, burnin=500, listw=listw_NY)
summary(res)
ecar1fw <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY, weights=POP8, family="CAR", method="eigen")
summary(ecar1fw)
res <- MCMCsamp(ecar1fw, mcmc=5000, burnin=500, listw=listw_NY)
summary(res)
esar0 <- errorsarlm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY)
summary(esar0)
res <- MCMCsamp(esar0, mcmc=5000, burnin=500, listw=listw_NY)
summary(res)
esar1 <- errorsarlm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY, etype="emixed")
summary(esar1)
res <- MCMCsamp(esar1, mcmc=5000, burnin=500, listw=listw_NY)
summary(res)
lsar0 <- lagsarlm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY)
summary(lsar0)
res <- MCMCsamp(lsar0, mcmc=5000, burnin=500, listw=listw_NY)
summary(res)
lsar1 <- lagsarlm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY, type="mixed")
summary(lsar1)
res <- MCMCsamp(lsar1, mcmc=5000, burnin=500, listw=listw_NY)
summary(res)
ssar0 <- sacsarlm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY)
summary(ssar0)
res <- MCMCsamp(ssar0, mcmc=5000, burnin=500, listw=listw_NY)
summary(res)
ssar1 <- sacsarlm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY, type="sacmixed")
summary(ssar1)
res <- MCMCsamp(ssar1, mcmc=5000, burnin=500, listw=listw_NY)
summary(res)
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{spatial}
