## ---------------------------------------------------------------------------------------
library(spdep)
NY8 <- as(sf::st_read(system.file("shapes/NY8_utm18.shp", package="spData")), "Spatial")
NY_nb <- read.gal(system.file("weights/NY_nb.gal", package="spData"), region.id=as.character(as.integer(row.names(NY8))-1L))

## ---------------------------------------------------------------------------------------
Syracuse <- NY8[NY8$AREANAME == "Syracuse city",]
Sy0_nb <- subset(NY_nb, NY8$AREANAME == "Syracuse city")
summary(Sy0_nb)

## ---------------------------------------------------------------------------------------
class(Syracuse)
Sy1_nb <- poly2nb(Syracuse)
isTRUE(all.equal(Sy0_nb, Sy1_nb, check.attributes=FALSE))

## ---------------------------------------------------------------------------------------
Sy2_nb <- poly2nb(Syracuse, queen=FALSE)
isTRUE(all.equal(Sy0_nb, Sy2_nb, check.attributes=FALSE))

## ----echo=FALSE,eval=TRUE---------------------------------------------------------------
run <- require("sp", quiet=TRUE)

## ----echo=run---------------------------------------------------------------------------
oopar <- par(mfrow=c(1,2), mar=c(3,3,1,1)+0.1)
plot(Syracuse, border="grey60")
plot(Sy0_nb, coordinates(Syracuse), add=TRUE, pch=19, cex=0.6)
text(bbox(Syracuse)[1,1], bbox(Syracuse)[2,2], labels="a)", cex=0.8)
plot(Syracuse, border="grey60")
plot(Sy0_nb, coordinates(Syracuse), add=TRUE, pch=19, cex=0.6)
plot(diffnb(Sy0_nb, Sy2_nb, verbose=FALSE), coordinates(Syracuse),
  add=TRUE, pch=".", cex=0.6, lwd=2)
text(bbox(Syracuse)[1,1], bbox(Syracuse)[2,2], labels="b)", cex=0.8)
par(oopar)

## ----eval=FALSE-------------------------------------------------------------------------
#  library(spgrass6)
#  writeVECT6(Syracuse, "SY0")
#  contig <- vect2neigh("SY0")
#  Sy3_nb <- sn2listw(contig)$neighbours
#  isTRUE(all.equal(Sy3_nb, Sy2_nb, check.attributes=FALSE))

## ----echo=run---------------------------------------------------------------------------
coords <- coordinates(Syracuse)
IDs <- row.names(as(Syracuse, "data.frame"))
#FIXME library(tripack)
Sy4_nb <- tri2nb(coords, row.names=IDs)
if (require(dbscan, quietly=TRUE)) {
  Sy5_nb <- graph2nb(soi.graph(Sy4_nb, coords), row.names=IDs)
} else Sy5_nb <- NULL
Sy6_nb <- graph2nb(gabrielneigh(coords), row.names=IDs)
Sy7_nb <- graph2nb(relativeneigh(coords), row.names=IDs)

## ----echo=run---------------------------------------------------------------------------
oopar <- par(mfrow=c(2,2), mar=c(1,1,1,1)+0.1)
plot(Syracuse, border="grey60")
plot(Sy4_nb, coords, add=TRUE, pch=".")
text(bbox(Syracuse)[1,1], bbox(Syracuse)[2,2], labels="a)", cex=0.8)
plot(Syracuse, border="grey60")
if (!is.null(Sy5_nb)) {
  plot(Sy5_nb, coords, add=TRUE, pch=".")
  text(bbox(Syracuse)[1,1], bbox(Syracuse)[2,2], labels="b)", cex=0.8)
}
plot(Syracuse, border="grey60")
plot(Sy6_nb, coords, add=TRUE, pch=".")
text(bbox(Syracuse)[1,1], bbox(Syracuse)[2,2], labels="c)", cex=0.8)
plot(Syracuse, border="grey60")
plot(Sy7_nb, coords, add=TRUE, pch=".")
text(bbox(Syracuse)[1,1], bbox(Syracuse)[2,2], labels="d)", cex=0.8)
par(oopar)

## ----echo=run---------------------------------------------------------------------------
nb_l <- list(Triangulation=Sy4_nb, Gabriel=Sy6_nb,
  Relative=Sy7_nb)
if (!is.null(Sy5_nb)) nb_l <- c(nb_l, list(SOI=Sy5_nb))
sapply(nb_l, function(x) is.symmetric.nb(x, verbose=FALSE, force=TRUE))
sapply(nb_l, function(x) n.comp.nb(x)$nc)

## ----echo=run---------------------------------------------------------------------------
Sy8_nb <- knn2nb(knearneigh(coords, k=1), row.names=IDs)
Sy9_nb <- knn2nb(knearneigh(coords, k=2), row.names=IDs)
Sy10_nb <- knn2nb(knearneigh(coords, k=4), row.names=IDs)
nb_l <- list(k1=Sy8_nb, k2=Sy9_nb, k4=Sy10_nb)
sapply(nb_l, function(x) is.symmetric.nb(x, verbose=FALSE, force=TRUE))
sapply(nb_l, function(x) n.comp.nb(x)$nc)

## ----echo=run---------------------------------------------------------------------------
oopar <- par(mfrow=c(1,3), mar=c(1,1,1,1)+0.1)
plot(Syracuse, border="grey60")
plot(Sy8_nb, coords, add=TRUE, pch=".")
text(bbox(Syracuse)[1,1], bbox(Syracuse)[2,2], labels="a)", cex=0.8)
plot(Syracuse, border="grey60")
plot(Sy9_nb, coords, add=TRUE, pch=".")
text(bbox(Syracuse)[1,1], bbox(Syracuse)[2,2], labels="b)", cex=0.8)
plot(Syracuse, border="grey60")
plot(Sy10_nb, coords, add=TRUE, pch=".")
text(bbox(Syracuse)[1,1], bbox(Syracuse)[2,2], labels="c)", cex=0.8)
par(oopar)

## ----echo=run---------------------------------------------------------------------------
dsts <- unlist(nbdists(Sy8_nb, coords))
summary(dsts)
max_1nn <- max(dsts)
max_1nn
Sy11_nb <- dnearneigh(coords, d1=0, d2=0.75*max_1nn, row.names=IDs)
Sy12_nb <- dnearneigh(coords, d1=0, d2=1*max_1nn, row.names=IDs)
Sy13_nb <- dnearneigh(coords, d1=0, d2=1.5*max_1nn, row.names=IDs)
nb_l <- list(d1=Sy11_nb, d2=Sy12_nb, d3=Sy13_nb)
sapply(nb_l, function(x) is.symmetric.nb(x, verbose=FALSE, force=TRUE))
sapply(nb_l, function(x) n.comp.nb(x)$nc)

## ----echo=run---------------------------------------------------------------------------
oopar <- par(mfrow=c(1,3), mar=c(1,1,1,1)+0.1)
plot(Syracuse, border="grey60")
plot(Sy11_nb, coords, add=TRUE, pch=".")
text(bbox(Syracuse)[1,1], bbox(Syracuse)[2,2], labels="a)", cex=0.8)
plot(Syracuse, border="grey60")
plot(Sy12_nb, coords, add=TRUE, pch=".")
text(bbox(Syracuse)[1,1], bbox(Syracuse)[2,2], labels="b)", cex=0.8)
plot(Syracuse, border="grey60")
plot(Sy13_nb, coords, add=TRUE, pch=".")
text(bbox(Syracuse)[1,1], bbox(Syracuse)[2,2], labels="c)", cex=0.8)
par(oopar)

## ----echo=run---------------------------------------------------------------------------
dS <- c(0.75, 1, 1.5)*max_1nn
res <- sapply(nb_l, function(x) table(card(x)))
mx <- max(card(Sy13_nb))
res1 <- matrix(0, ncol=(mx+1), nrow=3)
rownames(res1) <- names(res)
colnames(res1) <- as.character(0:mx)
res1[1, names(res$d1)] <- res$d1
res1[2, names(res$d2)] <- res$d2
res1[3, names(res$d3)] <- res$d3
library(RColorBrewer)
pal <- grey.colors(3, 0.95, 0.55, 2.2)
# RSB quietening greys
barplot(res1, col=pal, beside=TRUE, legend.text=FALSE, xlab="numbers of neighbours", ylab="tracts")
legend("topright", legend=format(dS, digits=1), fill=pal, bty="n", cex=0.8, title="max. distance")

## ----echo=run---------------------------------------------------------------------------
dsts0 <- unlist(nbdists(NY_nb, coordinates(NY8)))
summary(dsts0)

## ----echo=run---------------------------------------------------------------------------
Sy0_nb_lags <- nblag(Sy0_nb, maxlag=9)
names(Sy0_nb_lags) <- c("first", "second", "third", "fourth", "fifth", "sixth", "seventh", "eighth", "ninth")
res <- sapply(Sy0_nb_lags, function(x) table(card(x)))
mx <- max(unlist(sapply(Sy0_nb_lags, function(x) card(x))))
nn <- length(Sy0_nb_lags)
res1 <- matrix(0, ncol=(mx+1), nrow=nn)
rownames(res1) <- names(res)
colnames(res1) <- as.character(0:mx)
for (i in 1:nn) res1[i, names(res[[i]])] <- res[[i]]
res1

## ----echo=run---------------------------------------------------------------------------
cell2nb(7, 7, type="rook", torus=TRUE)
cell2nb(7, 7, type="rook", torus=FALSE)

## ----echo=run---------------------------------------------------------------------------
data(meuse.grid)
coordinates(meuse.grid) <- c("x", "y")
gridded(meuse.grid) <- TRUE
dst <- max(slot(slot(meuse.grid, "grid"), "cellsize"))
mg_nb <- dnearneigh(coordinates(meuse.grid), 0, dst)
mg_nb
table(card(mg_nb))

