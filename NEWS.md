# Version 1.2-3 (development)

* replace deprecated S-compatibility macros `DOUBLE_`

* #81 improved `dnearneigh()` help page.

* #79 remove `"htest"` class from `LOSH.mc()` output object.

* Added GA SI article to citations.

# Version 1.2-2 (2022-01-28)

* Replace `rainbow()` by `hcl.colors(..., "Set 2")` in `plot.skater()`.

* Add link to R-sig-geo thread on `EBlocal()` NaN estimates when many counts are zero on help page.

* Revise and add documentation for object returned by `localC_perm()` #68 #72 #73 #74 #75 #76

* `localmoran.sad()`, `localmoran.exact()` and `localmoran.exact.alt()` will now use multiple compute nodes if needed; if `Omega` is used, multiple cores may need more memory #77

* For **s2** > 1.0.7, use indexed distances in `dnearneigh()` https://github.com/r-spatial/s2/pull/162.

# Version 1.2-1 (2022-01-05)

* Switching deprecated functions moved to **spatialreg** to defunct.

# Version 1.1-13 (2021-12-14)

* Recent changes in `poly2nb()` had reduced and most recently (1.1-8) removed the use of `snap=` in finding candidate neighbours; many thanks to Matilda Brown for a clear and well-documented issue #65 

* Add local Geary's C #66 thanks to Josiah Parry, discussion on further work on #68

* `localmoran_perm()` returns both look-up and folded rank p-values

# Version 1.1-12 (2021-11-09)

* In `poly2nb()`, reverted removal of legacy interpreted overlapping envelope code for sp objects that cannot be coerced to sf without **rgeos**.

* Add Fortran character handling `USE_FC_LEN_T` WRE §6.6.1.

* Checks OK with forthcoming **deldir** 1.0-0.

* Fixes #62 clarifying `dnearneigh()` help page

# Version 1.1-11 (2021-09-07)

* `knearneigh()` and `nbdists()`; added prototype adaptation to **s2** for unprojected coordinates, used if `sf_use_s2()` is `TRUE` which became the default for **sf** 1.0.0 https://github.com/r-spatial/s2/issues/125. These are activated by default.

* `dnearneigh()` can choose the prototype **s2** approach if `sf_use_s2()` is `TRUE` and `use_s2=TRUE` for unprojected coordinates; from https://github.com/r-spatial/s2/issues/125 it seems that distance thresholds at present use brute-force rather than spatial indexing. Use is not activated by default.

* `poly2nb()` now uses `sf::st_intersects()` to find candidate neighbours unless `findInBounds=` is not NULL. With spatial indexing, this is very fast and scales well for large data sets. If `sf_use_s2()` is `TRUE`, `sf::st_intersects()` passes the geometries to `s2::s2_intersects_matrix()`, which also uses spatial indexing and is very fast, scaling well for large data sets.

* `localmoran()` and `localmoran_perm()` return cluster quadrants in an attribute for three splits, zeros, means and medians on the variable of interest and its spatial lag.

* `localmoran_perm()` returns the skewness and kurtosis of the permutation samples.


# Version 1.1-8 (2021-05-23)

* #55 related to #20 and cycling order in setting up grids provoked re-design of interface to `cell2nb()`, with passing of `"GridTopology"` or `"SpatialGrid"` objects as unnamed first or `x=` argument. Coerce `"RasterLayer"` or similar **raster**, **terra** or **stars** objects to **sp** class objects first if need be.

* In working with renewing the arguments to `cell2nb()`, it was useful to add **tinytest** support, which is now present for this function and may be extended to other functions for creating `"nb"` objects.

* #58 contributed by Jeff Sauer and Levi Wolf (from https://doi.org/10.31219/osf.io/ugkhp) providing conditional standard deviates for local Moran's I

* Error in assignment to matrix detected by CRAN check in SIDS vignette, section on median polish

# Version 1.1-7 (2021-04-03)

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
