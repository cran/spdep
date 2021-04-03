## ---- eval=FALSE------------------------------------------------------------------------
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

## ---- echo=FALSE------------------------------------------------------------------------
rv <- R.Version()
dothis <- FALSE
if (rv$major > "3" || (rv$major == "3" && !(rv$minor >= "3.0"))) dothis=TRUE

## ---- echo=dothis, eval=dothis----------------------------------------------------------
if (!suppressPackageStartupMessages(require(sf, quietly=TRUE))) {
  message("install the sf package")
  dothis <- FALSE
}
if (dothis) sf_extSoftVersion()

## ---- echo=dothis, eval=dothis----------------------------------------------------------
NY8_sf <- st_read(system.file("shapes/NY8_bna_utm18.gpkg", package="spData"), quiet=TRUE)
table(st_is_valid(NY8_sf))

## ---- echo=dothis, eval=dothis----------------------------------------------------------
suppressPackageStartupMessages(library(spdep))
reps <- 10
eps <- sqrt(.Machine$double.eps)
system.time(for(i in 1:reps) NY8_sf_1_nb <- poly2nb(NY8_sf, queen=TRUE, snap=eps))/reps

## ---- echo=dothis, eval=dothis----------------------------------------------------------
system.time(for(i in 1:reps) NY8_sf_1_nb_STR <- poly2nb(NY8_sf, queen=TRUE, snap=eps, small_n=250))/reps

## ---- echo=dothis, eval=dothis----------------------------------------------------------
NY8_sf_1_nb

## ---- echo=dothis, eval=dothis----------------------------------------------------------
all.equal(NY8_sf_1_nb, NY8_sf_1_nb_STR, check.attributes=FALSE)

## ---- echo=dothis, eval=dothis----------------------------------------------------------
NY8_sf_old <- st_read(system.file("shapes/NY8_utm18.shp", package="spData"), quiet=TRUE)
table(st_is_valid(NY8_sf_old))

## ---- echo=dothis, eval=dothis----------------------------------------------------------
try(NY8_sf_old_1_nb <- poly2nb(NY8_sf_old), silent = TRUE)
all.equal(NY8_sf_old_1_nb, NY8_sf_1_nb, check.attributes=FALSE)

## ---- echo=dothis, eval=dothis----------------------------------------------------------
NY8_sf_old_buf <- st_buffer(NY8_sf_old, dist=0)
table(st_is_valid(NY8_sf_old_buf))

## ---- echo=dothis, eval=dothis----------------------------------------------------------
try(NY8_sf_old_1_nb_buf <- poly2nb(NY8_sf_old_buf), silent = TRUE)
all.equal(NY8_sf_old_1_nb_buf, NY8_sf_1_nb, check.attributes=FALSE)

## ---- echo=dothis, eval=dothis----------------------------------------------------------
all.equal(NY8_sf_old_1_nb_buf, NY8_sf_old_1_nb, check.attributes=FALSE)

## ---- echo=dothis, eval=dothis----------------------------------------------------------
NY8_ct_sf <- st_centroid(st_geometry(NY8_sf), of_largest_polygon=TRUE)

## ---- echo=dothis, eval=dothis----------------------------------------------------------
NY8_pos_sf <- st_point_on_surface(st_geometry(NY8_sf))

## ---- echo=dothis, eval=dothis----------------------------------------------------------
if (unname(sf_extSoftVersion()["GEOS"] >= "3.9.0")) 
    NY8_cic_sf <- st_cast(st_inscribed_circle(st_geometry(NY8_sf), nQuadSegs=0), "POINT")[(1:(2*nrow(NY8_sf)) %% 2) != 0]

## ---- echo=dothis, eval=dothis----------------------------------------------------------
st_is_longlat(NY8_ct_sf)

## ---- echo=dothis, eval=dothis----------------------------------------------------------
suppressPackageStartupMessages(require(deldir))
NY84_nb <- tri2nb(NY8_ct_sf)
if (require(dbscan, quietly=TRUE)) {
  NY85_nb <- graph2nb(soi.graph(NY84_nb, NY8_ct_sf))
} else NY85_nb <- NULL
NY86_nb <- graph2nb(gabrielneigh(NY8_ct_sf))
NY87_nb <- graph2nb(relativeneigh(NY8_ct_sf))

## ---- echo=dothis, eval=dothis----------------------------------------------------------
system.time(for (i in 1:reps) NY88_nb_sf <- knn2nb(knearneigh(NY8_ct_sf, k=1)))/reps

## ---- echo=dothis, eval=dothis----------------------------------------------------------
system.time(for (i in 1:reps) NY89_nb_sf <- knn2nb(knearneigh(NY8_ct_sf, k=1, use_kd_tree=FALSE)))/reps

## ---- echo=dothis, eval=dothis----------------------------------------------------------
dsts <- unlist(nbdists(NY88_nb_sf, NY8_ct_sf))
summary(dsts)
max_1nn <- max(dsts)

## ---- echo=dothis, eval=dothis----------------------------------------------------------
system.time(for (i in 1:reps) NY810_nb <- dnearneigh(NY8_ct_sf, d1=0, d2=0.75*max_1nn))/reps

## ---- echo=dothis, eval=dothis----------------------------------------------------------
system.time(for (i in 1:reps) NY811_nb <- dnearneigh(NY8_ct_sf, d1=0, d2=0.75*max_1nn, use_kd_tree=FALSE))/reps

