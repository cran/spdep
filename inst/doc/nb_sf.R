## ----eval=FALSE-------------------------------------------------------------------------
#  library(sf)
#  sf_bna <- st_read("t8_36.bna", stringsAsFactors=FALSE)
#  table(st_is_valid(sf_bna))
#  sf_bna$AREAKEY <- gsub("\\.", "", sf_bna$Primary.ID)
#  data(NY_data, package="spData")
#  key <- as.character(nydata$AREAKEY)
#  sf_bna1 <- sf_bna[match(key, sf_bna$AREAKEY), c("AREAKEY")]
#  sf_bna2 <- merge(sf_bna1, nydata, by="AREAKEY")
#  sf_bna2_utm18 <- st_transform(sf_bna2, "+proj=utm +zone=18 +datum=NAD27")
#  st_write(sf_bna2_utm18, "NY8_bna_utm18.gpkg")

## ----echo=FALSE-------------------------------------------------------------------------
rv <- R.Version()
dothis <- FALSE
if (rv$major > "3" || (rv$major == "3" && !(rv$minor >= "3.0"))) dothis=TRUE

## ----echo=dothis, eval=dothis-----------------------------------------------------------
if (!suppressPackageStartupMessages(require(sf, quietly=TRUE))) {
  message("install the sf package")
  dothis <- FALSE
}
if (dothis) sf_extSoftVersion()

## ----echo=dothis, eval=dothis-----------------------------------------------------------
NY8_sf <- st_read(system.file("shapes/NY8_bna_utm18.gpkg", package="spData"), quiet=TRUE)
table(st_is_valid(NY8_sf))

## ----echo=dothis, eval=dothis-----------------------------------------------------------
suppressPackageStartupMessages(library(spdep))
reps <- 10
eps <- sqrt(.Machine$double.eps)
system.time(for(i in 1:reps) NY8_sf_1_nb <- poly2nb(NY8_sf, queen=TRUE, snap=eps))/reps

## ----echo=dothis, eval=dothis-----------------------------------------------------------
NY8_sf_1_nb

## ----echo=dothis, eval=dothis-----------------------------------------------------------
NY8_sf_old <- st_read(system.file("shapes/NY8_utm18.shp", package="spData"), quiet=TRUE)
table(st_is_valid(NY8_sf_old))

## ----echo=dothis, eval=dothis-----------------------------------------------------------
try(NY8_sf_old_1_nb <- poly2nb(NY8_sf_old), silent = TRUE)
all.equal(NY8_sf_old_1_nb, NY8_sf_1_nb, check.attributes=FALSE)

## ----echo=dothis, eval=dothis-----------------------------------------------------------
NY8_sf_old_val <- st_make_valid(NY8_sf_old, dist=0)
table(st_is_valid(NY8_sf_old_val))

## ----echo=dothis, eval=dothis-----------------------------------------------------------
class(st_geometry(NY8_sf_old))

## ----echo=dothis, eval=dothis-----------------------------------------------------------
class(st_geometry(NY8_sf_old_val))

## ----echo=dothis, eval=dothis-----------------------------------------------------------
table(sapply(st_geometry(NY8_sf_old_val), function(x) class(x)[[2]]))

## ----echo=dothis, eval=dothis-----------------------------------------------------------
NY8_sf_old_val <- st_collection_extract(NY8_sf_old_val, "POLYGON")
table(sapply(st_geometry(NY8_sf_old_val), function(x) class(x)[[2]]))

## ----echo=dothis, eval=dothis-----------------------------------------------------------
try(NY8_sf_old_1_nb_val <- poly2nb(NY8_sf_old_val), silent = TRUE)
all.equal(NY8_sf_old_1_nb_val, NY8_sf_1_nb, check.attributes=FALSE)

## ----echo=dothis, eval=dothis-----------------------------------------------------------
all.equal(NY8_sf_old_1_nb_val, NY8_sf_old_1_nb, check.attributes=FALSE)

## ----echo=dothis, eval=dothis-----------------------------------------------------------
NY8_ct_sf <- st_centroid(st_geometry(NY8_sf), of_largest_polygon=TRUE)

## ----echo=dothis, eval=dothis-----------------------------------------------------------
NY8_pos_sf <- st_point_on_surface(st_geometry(NY8_sf))

## ----echo=dothis, eval=dothis-----------------------------------------------------------
if (unname(sf_extSoftVersion()["GEOS"] >= "3.9.0")) 
    NY8_cic_sf <- st_cast(st_inscribed_circle(st_geometry(NY8_sf), nQuadSegs=0), "POINT")[(1:(2*nrow(NY8_sf)) %% 2) != 0]

## ----echo=dothis, eval=dothis-----------------------------------------------------------
st_is_longlat(NY8_ct_sf)

## ----echo=dothis, eval=dothis-----------------------------------------------------------
suppressPackageStartupMessages(require(deldir))
NY84_nb <- tri2nb(NY8_ct_sf)
if (require(dbscan, quietly=TRUE)) {
  NY85_nb <- graph2nb(soi.graph(NY84_nb, NY8_ct_sf))
} else NY85_nb <- NULL
NY86_nb <- graph2nb(gabrielneigh(NY8_ct_sf))
NY87_nb <- graph2nb(relativeneigh(NY8_ct_sf))

## ----echo=dothis, eval=dothis-----------------------------------------------------------
system.time(for (i in 1:reps) NY88_nb_sf <- knn2nb(knearneigh(NY8_ct_sf, k=1)))/reps

## ----echo=dothis, eval=dothis-----------------------------------------------------------
system.time(for (i in 1:reps) NY89_nb_sf <- knn2nb(knearneigh(NY8_ct_sf, k=1, use_kd_tree=FALSE)))/reps

## ----echo=dothis, eval=dothis-----------------------------------------------------------
dsts <- unlist(nbdists(NY88_nb_sf, NY8_ct_sf))
summary(dsts)
max_1nn <- max(dsts)

## ----echo=dothis, eval=dothis-----------------------------------------------------------
system.time(for (i in 1:reps) NY810_nb <- dnearneigh(NY8_ct_sf, d1=0, d2=0.75*max_1nn))/reps

## ----echo=dothis, eval=dothis-----------------------------------------------------------
system.time(for (i in 1:reps) NY811_nb <- dnearneigh(NY8_ct_sf, d1=0, d2=0.75*max_1nn, use_kd_tree=FALSE))/reps

## ----echo=dothis, eval=dothis-----------------------------------------------------------
pts_ll <- st_transform(NY8_ct_sf, "OGC:CRS84")
st_is_longlat(pts_ll)

## ----echo=dothis, eval=dothis-----------------------------------------------------------
(old_use_s2 <- sf_use_s2())

## ----echo=dothis, eval=dothis-----------------------------------------------------------
sf_use_s2(TRUE)
system.time(for (i in 1:reps) pts_ll1_nb <- knn2nb(knearneigh(pts_ll, k=6)))/reps

## ----echo=dothis, eval=dothis-----------------------------------------------------------
sf_use_s2(FALSE)
system.time(for (i in 1:reps) pts_ll2_nb <- knn2nb(knearneigh(pts_ll, k=6)))/reps

## ----echo=dothis, eval=dothis-----------------------------------------------------------
all.equal(pts_ll1_nb, pts_ll2_nb, check.attributes=FALSE)

## ----echo=dothis, eval=dothis-----------------------------------------------------------
pts_ll1_nb[[52]]
pts_ll2_nb[[52]]
pts_ll1_nb[[124]]
pts_ll2_nb[[124]]

## ----echo=dothis, eval=dothis-----------------------------------------------------------
sf_use_s2(old_use_s2)

## ----echo=dothis, eval=dothis-----------------------------------------------------------
max_1nn_ll <- max(unlist(nbdists(knn2nb(knearneigh(pts_ll, k=1)), pts_ll)))

## ----echo=dothis, eval=dothis-----------------------------------------------------------
args(dnearneigh)

## ----echo=dothis, eval=dothis-----------------------------------------------------------
if (packageVersion("s2") > "1.0.7") {
  system.time(for (i in 1:(reps/5)) pts_ll3_nb <- dnearneigh(pts_ll, d1=0,
      d2=0.75*max_1nn_ll))/(reps/5)
}

## ----echo=dothis, eval=dothis-----------------------------------------------------------
system.time(for (i in 1:(reps/5)) pts_ll5_nb <- dnearneigh(pts_ll, d1=0, d2=0.75*max_1nn_ll, dwithin=FALSE))/(reps/5)

## ----echo=dothis, eval=dothis-----------------------------------------------------------
if (packageVersion("s2") > "1.0.7") all.equal(pts_ll3_nb, pts_ll5_nb, check.attributes=FALSE)

## ----echo=dothis, eval=dothis-----------------------------------------------------------
if (packageVersion("s2") > "1.0.7") {
  system.time(for (i in 1:(reps/5)) pts_ll3a_nb <- dnearneigh(pts_ll, d1=5,
      d2=0.75*max_1nn_ll, dwithin=FALSE))/(reps/5)
}

## ----echo=dothis, eval=dothis-----------------------------------------------------------
if (packageVersion("s2") > "1.0.7") {
    system.time(for (i in 1:(reps/5)) pts_ll5a_nb <- dnearneigh(pts_ll, d1=5,
        d2=0.75*max_1nn_ll))/(reps/5)
}

## ----echo=dothis, eval=dothis-----------------------------------------------------------
if (packageVersion("s2") > "1.0.7") all.equal(pts_ll3a_nb, pts_ll5a_nb, check.attributes=FALSE)

## ----echo=dothis, eval=dothis-----------------------------------------------------------
system.time(for (i in 1:reps) pts_ll6_nb <- dnearneigh(pts_ll, d1=0, d2=0.75*max_1nn_ll, use_s2=FALSE))/reps

## ----echo=dothis, eval=dothis-----------------------------------------------------------
all.equal(pts_ll5_nb, pts_ll6_nb, check.attributes=FALSE)

## ----echo=dothis, eval=dothis-----------------------------------------------------------
system.time(for (i in 1:reps) pts_ll6a_nb <- dnearneigh(pts_ll, d1=5, d2=0.75*max_1nn_ll, use_s2=FALSE))/reps

## ----echo=dothis, eval=dothis-----------------------------------------------------------
if (packageVersion("s2") > "1.0.7") all.equal(pts_ll5a_nb, pts_ll6a_nb, check.attributes=FALSE)

## ----echo=dothis, eval=dothis-----------------------------------------------------------
NY8_sf_ll <- st_transform(NY8_sf, "OGC:CRS84")
st_is_longlat(NY8_sf_ll)

## ----echo=dothis, eval=dothis-----------------------------------------------------------
system.time(for (i in 1:reps) NY8_sf_1_nb_ll <- poly2nb(NY8_sf_ll, queen=TRUE, snap=eps))/reps

## ----echo=dothis, eval=dothis-----------------------------------------------------------
all.equal(NY8_sf_1_nb, NY8_sf_1_nb_ll, check.attributes=FALSE)

