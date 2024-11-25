## ----setup, include=FALSE---------------------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---------------------------------------------------------------------------------------
library(spdep)
args(moran.test)

## ---------------------------------------------------------------------------------------
eigen(0)$values

## ---------------------------------------------------------------------------------------
(GDAL37 <- as.numeric_version(unname(sf_extSoftVersion()["GDAL"])) >= "3.7.0")

## ---------------------------------------------------------------------------------------
file <- "etc/shapes/GB_2024_Wales_50m.gpkg.zip"
zipfile <- system.file(file, package="spdep")
if (GDAL37) {
    w50m <- st_read(zipfile)
} else {
    td <- tempdir()
    bn <- sub(".zip", "", basename(file), fixed=TRUE)
    target <- unzip(zipfile, files=bn, exdir=td)
    w50m <- st_read(target)
}

## ---------------------------------------------------------------------------------------
(w50m |> poly2nb(row.names=as.character(w50m$Constituency)) -> nb_W_50m)

## ---------------------------------------------------------------------------------------
attr(nb_W_50m, "ncomp")$comp.id |>table() |> table()

## ---------------------------------------------------------------------------------------
ynys_mon <- w50m$Constituency == "Ynys Môn"
pts <- st_point_on_surface(st_geometry(w50m))
opar <- par(mfrow=c(1, 2))
plot(st_geometry(w50m), border="grey75")
plot(st_geometry(w50m)[ynys_mon], add=TRUE)
plot(st_geometry(w50m)[card(nb_W_50m) == 0L], add=TRUE, border="transparent", col="wheat1")
plot(st_geometry(w50m), border="grey75")
plot(nb_W_50m, pts, add=TRUE)
par(opar)

## ---------------------------------------------------------------------------------------
dym <- c(st_distance(w50m[ynys_mon,], w50m))
sort(dym)[1:12]

## ---------------------------------------------------------------------------------------
(nb_W_50m_snap <- poly2nb(w50m, row.names=as.character(w50m$Constituency), snap=280))

## ---------------------------------------------------------------------------------------
plot(st_geometry(w50m), border="grey75")
plot(nb_W_50m_snap, pts, add=TRUE)

## ---------------------------------------------------------------------------------------
attr(nb_W_50m_snap, "region.id")[nb_W_50m_snap[[which(ynys_mon)]]]

## ---------------------------------------------------------------------------------------
(meet_criterion <- sum(dym <= units::set_units(280, "m")))

## ---------------------------------------------------------------------------------------
(cands <- attr(nb_W_50m, "region.id")[order(dym)[1:meet_criterion]])

## ---------------------------------------------------------------------------------------
(nb_W_50m_add <- addlinks1(nb_W_50m, from = cands[1], to = cands[2:meet_criterion]))

## ---------------------------------------------------------------------------------------
all.equal(nb_W_50m_add, nb_W_50m_snap, check.attributes=FALSE)

## ---------------------------------------------------------------------------------------
k2 <- knn2nb(knearneigh(pts, k=2), row.names=as.character(w50m$Constituency), sym=TRUE)
attr(k2, "region.id")[k2[[which(ynys_mon)]]]

## ---------------------------------------------------------------------------------------
(k6 <- knn2nb(knearneigh(pts, k=6), row.names=as.character(w50m$Constituency), sym=TRUE))

## ---------------------------------------------------------------------------------------
plot(st_geometry(w50m), border="grey75")
plot(k6, pts, add=TRUE)

## ---------------------------------------------------------------------------------------
o <- order(attr(k6, "ncomp")$comp.id)
image(t(nb2mat(k6, style="B")[o, rev(o)]), axes=FALSE, asp=1)

## ---------------------------------------------------------------------------------------
(k6a <- knn2nb(knearneigh(pts, k=6), row.names=as.character(w50m$Constituency)))

## ---------------------------------------------------------------------------------------
file <- "etc/shapes/GB_2024_southcoast_50m.gpkg.zip"
zipfile <- system.file(file, package="spdep")
if (GDAL37) {
    sc50m <- st_read(zipfile)
} else {
    td <- tempdir()
    bn <- sub(".zip", "", basename(file), fixed=TRUE)
    target <- unzip(zipfile, files=bn, exdir=td)
    sc50m <- st_read(target)
}

## ---------------------------------------------------------------------------------------
(nb_sc_50m <- poly2nb(sc50m, row.names=as.character(sc50m$Constituency)))

## ---------------------------------------------------------------------------------------
nc <- attr(nb_sc_50m, "ncomp")$comp.id
table(nc)

## ---------------------------------------------------------------------------------------
(sub2 <- attr(nb_sc_50m, "region.id")[nc == 2L])

## ---------------------------------------------------------------------------------------
pts <- st_point_on_surface(st_geometry(sc50m))
plot(st_geometry(sc50m), border="grey75")
plot(st_geometry(sc50m)[nc == 2L], border="orange", lwd=2, add=TRUE)
plot(nb_sc_50m, pts, add=TRUE)

## ---------------------------------------------------------------------------------------
1/range(eigen(cbind(c(0, 1), c(1, 0)))$values)
1/range(eigen(nb2mat(subset(nb_sc_50m, nc == 2L), style="W"))$values)

## ---------------------------------------------------------------------------------------
1/range(eigen(nb2mat(nb_sc_50m, style="W"))$values)

## ---------------------------------------------------------------------------------------
1/range(eigen(nb2mat(subset(nb_sc_50m, nc == 1L), style="W"))$values)

## ---------------------------------------------------------------------------------------
iowe <- match(sub2[1], attr(nb_sc_50m, "region.id"))
diowe <- c(st_distance(sc50m[iowe,], sc50m))
sort(diowe)[1:12]

## ---------------------------------------------------------------------------------------
ioww <- match(sub2[2], attr(nb_sc_50m, "region.id"))
dioww <- c(st_distance(sc50m[ioww,], sc50m))
sort(dioww)[1:12]

## ---------------------------------------------------------------------------------------
(meet_criterion <- sum(diowe <= units::set_units(5000, "m")))

## ---------------------------------------------------------------------------------------
(cands <- attr(nb_sc_50m, "region.id")[order(diowe)[1:meet_criterion]])

## ---------------------------------------------------------------------------------------
(nb_sc_50m_iowe <- addlinks1(nb_sc_50m, from = cands[1], to = cands[3:meet_criterion]))

## ---------------------------------------------------------------------------------------
(meet_criterion <- sum(dioww <= units::set_units(5000, "m")))

## ---------------------------------------------------------------------------------------
(cands <- attr(nb_sc_50m, "region.id")[order(dioww)[1:meet_criterion]])

## ---------------------------------------------------------------------------------------
(nb_sc_50m_iow <- addlinks1(nb_sc_50m_iowe, from = cands[2], to = cands[3:meet_criterion]))

## ---------------------------------------------------------------------------------------
pts <- st_point_on_surface(st_geometry(sc50m))
plot(st_geometry(sc50m), border="grey75")
plot(st_geometry(sc50m)[nc == 2L], border="orange", lwd=2, add=TRUE)
plot(nb_sc_50m_iow, pts, add=TRUE)

## ---------------------------------------------------------------------------------------
get.ZeroPolicyOption()

## ---------------------------------------------------------------------------------------
try(nb2listw(nb_W_50m))

## ---------------------------------------------------------------------------------------
set.ZeroPolicyOption(TRUE)

## ---------------------------------------------------------------------------------------
get.ZeroPolicyOption()

## ----eval=FALSE, echo=TRUE--------------------------------------------------------------
# (lw <- nb2listw(nb_W_50m))

## ----echo=FALSE-------------------------------------------------------------------------
# repeated to overcome CMD build latency
(lw <- nb2listw(nb_W_50m, zero.policy=get.ZeroPolicyOption()))

## ---------------------------------------------------------------------------------------
attr(lw, "zero.policy")

## ---------------------------------------------------------------------------------------
set.ZeroPolicyOption(FALSE)

## ---------------------------------------------------------------------------------------
get.NoNeighbourOption()
get.SubgraphOption()
get.SubgraphCeiling()

## ---------------------------------------------------------------------------------------
set.NoNeighbourOption(FALSE)
(w50m |> poly2nb(row.names=as.character(w50m$Constituency)) -> nb_W_50mz)

## ---------------------------------------------------------------------------------------
set.SubgraphOption(FALSE)
(w50m |> poly2nb(row.names=as.character(w50m$Constituency)) -> nb_W_50my)

## ---------------------------------------------------------------------------------------
str(attr(nb_W_50my, "ncomp"))

## ---------------------------------------------------------------------------------------
set.SubgraphOption(TRUE)
set.SubgraphCeiling(100L)
(w50m |> poly2nb(row.names=as.character(w50m$Constituency)) -> nb_W_50mx)

## ---------------------------------------------------------------------------------------
str(attr(nb_W_50mx, "ncomp"))

## ---------------------------------------------------------------------------------------
set.SubgraphCeiling(100000L)
set.NoNeighbourOption(TRUE)

## ---------------------------------------------------------------------------------------
file <- "etc/shapes/tokyo.gpkg.zip"
zipfile <- system.file(file, package="spdep")
if (GDAL37) {
    tokyo <- st_read(zipfile)
} else {
    td <- tempdir()
    bn <- sub(".zip", "", basename(file), fixed=TRUE)
    target <- unzip(zipfile, files=bn, exdir=td)
    tokyo <- st_read(target)
}

## ---------------------------------------------------------------------------------------
all(st_is_valid(tokyo))
tokyo <- st_make_valid(tokyo)

## ---------------------------------------------------------------------------------------
(nb_t0 <- poly2nb(tokyo, snap=sqrt(.Machine$double.eps)))

## ---------------------------------------------------------------------------------------
units::set_units(units::set_units(attr(nb_t0, "snap"), "m"), "nm")

## ---------------------------------------------------------------------------------------
(nb_t1 <- poly2nb(tokyo, snap=0.002))

## ---------------------------------------------------------------------------------------
units::set_units(units::set_units(attr(nb_t1, "snap"), "m"), "mm")

## ---------------------------------------------------------------------------------------
(nb_t2 <- poly2nb(tokyo))

## ---------------------------------------------------------------------------------------
units::set_units(units::set_units(attr(nb_t2, "snap"), "m"), "mm")

## ---------------------------------------------------------------------------------------
(nb_t3 <- poly2nb(st_transform(tokyo, "OGC:CRS84")))

## ---------------------------------------------------------------------------------------
attr(nb_t3, "snap")

## ---------------------------------------------------------------------------------------
(180 * 0.01) / (pi * 6378137)

