# Version 1.1-7 (development)

* Changes to continuous integration and vignettes.

* Error in `poly2nb(, queen=FALSE)` in **sf** grids (double counting of closed polygon start/end points), https://github.com/r-spatial/spdep/issues/50, thanks to Christopher Kenny.

* Adding local Moran and local G conditional permutation: `localmoran_perm()` and `localG_perm()`.

* Adding `nb2listwdist()` contributed by René Westerholt.

* Adding use of **sf** through GEOS to find polygon contiguity candidates in `poly2nb()` if geometry count >= 500 - uses intersections in polygon envelopes.

* #38, #53 removing **RANN**, adding **dbscan** suggestions for fast `dnearneigh()` and `knearneigh()` via `use_kd_tree=` argument for fast planar neighbour set finding in 2D and 3D. Affects `soi.graph()` too, which had used **RANN**.

* #54 avoid partial matching in `glist=` handling.

* Disambiguating **spdep** and **spatialreg** model output object class names prior to making **spdep** model fitting functions defunct.

# Version 1.1-5 (2020-06-29)

* Replacing broken geoda URLs, moving knitr to rmarkdown, work-around missing weights files in spData.


# Version 1.1-3 (2019-09-18)

* A small maintenance update to accommodate a forthcoming change in spData (a dataset used in an example in spdep from spData is changing its name;   the name had involved putting "x", "y" and "xyz" in the global environment through lazy loading a dataset).


# Version 1.1-2 (2019-04-05)

* Follow-up version of spdep with all the functions and 
  methods transferred to the spatialreg package marked 
  as deprecated, but still exported from spdep. Reverse 
  dependencies passing with the released version still 
  pass for me with this version.
