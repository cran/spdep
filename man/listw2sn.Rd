% Copyright 2001-2 by Roger S. Bivand
\name{listw2sn}
\alias{listw2sn}
\alias{spwdet}
\alias{log.spwdet}
\title{Spatial neighbour sparse representation}
\description{
The function makes a \code{spatial neighbour} object representation
(similar to the S-PLUS spatial statististics module representation of a
\code{listw} spatial weights object. The object is used in calculating the
Jacobian of the likelihood functions of spatial autoregressive models when
the "sparse" method is chosen in \code{errorsarlm} and \code{lagsarlm}.
}
\usage{
listw2sn(listw)
spwdet(sparseweights, rho)
log.spwdet(sparseweights, rho)

}
\arguments{
  \item{listw}{a \code{listw} object from for example \code{nb2listw}}
  \item{sparseweights}{a \code{spatial.neighbour} object}
  \item{rho}{spatial autoregressive parameter value}
}

\details{The methods used for finding the determinant directly by sparse
matrix techniques are given by R. Kelley Pace and R. P. Barry (1998),
and are implemented here using the sparse library by Kenneth Kundert
and Alberto Sangiovanni-Vincentelli; the same reference is used by the
S-PLUS spatial statistics module.}

\value{
A data frame with three columns, and with class \code{spatial.neighbour}:
  \item{from}{region number id for the start of the link (S-PLUS row.id)}
  \item{to}{region number id for the end of the link (S-PLUS col.id)}
  \item{weights}{weight for this link}
\code{log.spwdet} returns log det(I - rho * W).
}

\references{R. Kelley Pace and R. P. Barry (1998) Quick computations
of spatially autoregressive estimators, Geographical Analysis 29 (3)
232-247;
\url{http://www.spatial-statistics.com/pace_manuscripts/ga_ms_dir/pdf/fin_geo_analysis.pdf}.

Files: src/spAllocate.c, src/spBuild.c, src/spConfig.h,
src/spDefs.h,
src/spFactor.c, src/spMatrix.h, src/spUtils.c are by Kenneth Kundert
and Alberto Sangiovanni-Vincentelli, University of California, Berkeley,
and are from:

\url{http://www.netlib.org/sparse/index.html}

(Kenneth Kundert, Sparse Matrix Techniques, in Circuit Analysis,
Simulation and Design, Albert Ruehli (Ed.), North-Holland, 1986)

They are copyright (c) 1985,86,87,88 by Kenneth S. Kundert and the
University of California, with the following licence:

Permission to use, copy, modify, and distribute this software and
its documentation for any purpose and without fee is hereby granted,
provided that the copyright notices appear in all copies and
supporting documentation and that the authors and the University of
California are properly credited.  The authors and the University of
California make no representations as to the suitability of this
software for any purpose.  It is provided "as is", without express
or implied warranty.}

\author{Roger Bivand \email{Roger.Bivand@nhh.no}}

\seealso{\code{\link{nb2listw}}, \code{\link{errorsarlm}},
 \code{\link{lagsarlm}}}

\examples{
data(columbus)
col.listw <- nb2listw(col.gal.nb)
col.listw$neighbours[[1]]
col.listw$weights[[1]]
col.sn <- listw2sn(col.listw)
col.sn[col.sn[,1] == 1,]
rho <- seq(-0.8, 0.9, 0.1)
for (i in rho) print(paste("rho:", i, "log(det(I - rho*W))",
  log.spwdet(col.sn, i)), quote=TRUE)
}
\keyword{spatial}
