## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(message = FALSE, warning = FALSE)

## ----echo=TRUE-----------------------------------------------------------
library(spdep)
if (require(rgdal, quietly=TRUE)) {
  NY8 <- readOGR(system.file("shapes/NY8_utm18.shp", package="spData"))
} else {
  require(maptools, quietly=TRUE)
  NY8 <- readShapeSpatial(system.file("shapes/NY8_utm18.shp", package="spData"))
}
NY_nb <- read.gal(system.file("weights/NY_nb.gal", package="spData"), region.id=row.names(NY8))

## ----echo=TRUE,eval=TRUE-------------------------------------------------
nySFE <- SpatialFiltering(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, nb=NY_nb, style="W", verbose=FALSE)
nylmSFE <- lm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME+fitted(nySFE), data=NY8)
summary(nylmSFE)
nylm <- lm(Z~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8)
anova(nylm, nylmSFE)

## ----echo=TRUE-----------------------------------------------------------
NYlistwW <- nb2listw(NY_nb, style = "W")
set.seed(111)

## ----echo=TRUE,eval=FALSE------------------------------------------------
#  nyME <- ME(Cases~PEXPOSURE+PCTAGE65P+PCTOWNHOME, data=NY8, offset=log(POP8), family="poisson", listw=NYlistwW, alpha=0.5)

## ----echo=FALSE, eval=FALSE----------------------------------------------
#  bsfn <- system.file("etc/backstore/nyME_res.RData", package="spdep")
#  load(bsfn)

## ----echo=TRUE, eval=FALSE-----------------------------------------------
#  nyME
#  NY8$eigen_24 <- fitted(nyME)[,1]
#  NY8$eigen_223 <- fitted(nyME)[,2]

## ----results='asis',echo=FALSE, eval=FALSE-------------------------------
#  .iwidth <- 6
#  .iheight <- 4
#  .ipointsize <- 10
#  library(RColorBrewer)
#  #gry <- brewer.pal(9, "Greys")[-1]
#  spplot(NY8, c("eigen_24", "eigen_223"), col.regions=grey.colors(6, 0.95, 0.55, 2.2), cuts=5)

## ----echo=TRUE,eval=FALSE------------------------------------------------
#  nyglmME <- glm(Cases~PEXPOSURE+PCTAGE65P+PCTOWNHOME+offset(log(POP8))+fitted(nyME), data=NY8, family="poisson")
#  summary(nyglmME)
#  nyGLMp <- glm(Cases~PEXPOSURE+PCTAGE65P+PCTOWNHOME+offset(log(POP8)), data=NY8,family="poisson")
#  anova(nyGLMp, nyglmME, test="Chisq")

