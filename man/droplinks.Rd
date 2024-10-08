% Copyright 2001 by Roger S. Bivand
\name{droplinks}
\alias{droplinks}
\alias{addlinks1}
\title{Drop  and add links in a neighbours list}
\description{
 \code{droplinks} drops links to and from or just to a region from a neighbours list. The example corresponds to Fingleton's Table 1, (1999) p. 6, for lattices 5 to 19. \code{addlinks1} adds links from a single region to specified regions.
}
\usage{
droplinks(nb, drop, sym=TRUE)
addlinks1(nb, from, to, sym=TRUE)
}

\arguments{
  \item{nb}{a neighbours list object of class \code{nb}}
  \item{drop}{either a logical vector the length of \code{nb}, or a character vector of named regions corresponding to \code{nb}'s region.id attribute, or an integer vector of region numbers}
  \item{sym}{TRUE for removal of both "row" and "column" links, FALSE for only "row" links; when adding links, inserts links to the from region from the to regions}
  \item{from}{single from region for adding links, either a character vector of length 1 of the named from region corresponding to \code{nb}'s region.id attribute, or an integer vector of length 1 holding a region number}
  \item{to}{to regions, either a character vector of named from regions corresponding to \code{nb}'s region.id attribute, or an integer vector of region numbers}
}
\value{
 The function returns an object of class \code{nb} with a list of integer vectors containing neighbour region number ids. 
}
\references{B. Fingleton (1999) Spurious spatial regression: some Monte Carlo results with a spatial unit root and spatial cointegration, Journal of Regional Science 39, pp. 1--19.}
\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{is.symmetric.nb}}}

\examples{
\donttest{
rho <- c(0.2, 0.5, 0.95, 0.999, 1.0)
ns <- c(5, 7, 9, 11, 13, 15, 17, 19)
mns <- matrix(0, nrow=length(ns), ncol=length(rho))
rownames(mns) <- ns
colnames(mns) <- rho
mxs <- matrix(0, nrow=length(ns), ncol=length(rho))
rownames(mxs) <- ns
colnames(mxs) <- rho
for (i in 1:length(ns)) {
  nblist <- cell2nb(ns[i], ns[i])
  nbdropped <- droplinks(nblist, ((ns[i]*ns[i])+1)/2, sym=FALSE)
  listw <- nb2listw(nbdropped, style="W", zero.policy=TRUE)
  wmat <- listw2mat(listw)
  for (j in 1:length(rho)) {
    mat <- diag(ns[i]*ns[i]) - rho[j] * wmat
    res <- diag(solve(t(mat) \%*\% mat))
    mns[i,j] <- mean(res)
    mxs[i,j] <- max(res)
  }
}
print(mns)
print(mxs)

}

}
\keyword{spatial}
