# Calculates the autocovariate to be used in autonormal, autopoisson 
# or autologistic regression. Three distance-weighting schemes are available
#
# z is the response variable
# xy is the matrix of coordinates
# nbs is "neighbourhood size", selected by user; default is 1
# type defines the weighting scheme: 
#	"one" gives equal weight to all data points in the neighbourhood; 
#	"inverse" (the default) weights by inverse Euclidean distance;
#	"inverse.squared" weights by the square of "inverse"
#
# by Carsten F. Dormann, 04.11.2005, carsten.dormann@ufz.de
# Re-implementation allowing list representation
# Roger Bivand 28.11.2005
#

autocov_dist <- function(z, xy, nbs=1, type="inverse", zero.policy=FALSE,
   style="W", lonlat=FALSE) {
   if (type=="one") expo <- 0
   if (type=="inverse") expo <- 1
   if (type=="inverse.squared") expo <- 2
   nb <- dnearneigh(xy, 0, nbs, lonlat=lonlat)
   if (any(card(nb) == 0)) warning(paste("With value", nbs,
      "some points have no neighbours"))
   nbd <- nbdists(nb, xy, lonlat=lonlat)
   if (expo == 0) lw <- nb2listw(nb, style=style, zero.policy=zero.policy)
   else {
      gl <- lapply(nbd, function(x) 1/(x^expo))
      lw <- nb2listw(nb, glist=gl, style=style, zero.policy=zero.policy)
   }
   lag(lw, z, zero.policy=zero.policy)
}


