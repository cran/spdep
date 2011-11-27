pkgname <- "spdep"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('spdep')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("EBImoran.mc")
### * EBImoran.mc

flush(stderr()); flush(stdout())

### Name: EBImoran.mc
### Title: Permutation test for empirical Bayes index
### Aliases: EBImoran.mc EBImoran
### Keywords: spatial

### ** Examples

example(nc.sids)
EBImoran.mc(nc.sids$SID74, nc.sids$BIR74,
 nb2listw(ncCC89_nb, style="B", zero.policy=TRUE), nsim=999, zero.policy=TRUE)
sids.p <- nc.sids$SID74 / nc.sids$BIR74
moran.mc(sids.p, nb2listw(ncCC89_nb, style="B", zero.policy=TRUE),
 nsim=999, zero.policy=TRUE)



cleanEx()
nameEx("EBest")
### * EBest

flush(stderr()); flush(stdout())

### Name: EBest
### Title: Global Empirical Bayes estimator
### Aliases: EBest
### Keywords: spatial

### ** Examples

example(auckland)
res <- EBest(auckland$M77_85, 9*auckland$Und5_81)
attr(res, "parameters")
cols <- grey(6:2/7)
brks <- c(-Inf,2,2.5,3,3.5,Inf)
plot(auckland, col=cols[findInterval(res$estmm*1000, brks, all.inside=TRUE)])
legend("bottomleft", fill=cols, legend=leglabs(brks), bty="n")
title(main="Global moment estimator of infant mortality per 1000 per year")
data(huddersfield)
res <- EBest(huddersfield$cases, huddersfield$total, family="binomial")
round(res[,1:2],4)*100



cleanEx()
nameEx("EBlocal")
### * EBlocal

flush(stderr()); flush(stdout())

### Name: EBlocal
### Title: Local Empirical Bayes estimator
### Aliases: EBlocal
### Keywords: spatial

### ** Examples

example(auckland)
res <- EBlocal(auckland$M77_85,  9*auckland$Und5_81, auckland.nb)
brks <- c(-Inf,2,2.5,3,3.5,Inf)
cols <- grey(6:2/7)
plot(auckland, col=cols[findInterval(res$est*1000, brks, all.inside=TRUE)])
legend("bottomleft", fill=cols, legend=leglabs(brks), bty="n")
title(main="Local moment estimator of infant mortality per 1000 per year")



cleanEx()
nameEx("GMerrorsar")
### * GMerrorsar

flush(stderr()); flush(stdout())

### Name: GMerrorsar
### Title: Spatial simultaneous autoregressive error model estimation by
###   GMM
### Aliases: GMerrorsar residuals.gmsar deviance.gmsar coef.gmsar
###   fitted.gmsar print.gmsar summary.gmsar print.summary.gmsar
###   GMargminImage
### Keywords: spatial

### ** Examples

data(oldcol)
COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="W"), method="eigen")
summary(COL.errW.eig, Hausman=TRUE)
COL.errW.GM <- GMerrorsar(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="W"), returnHcov=TRUE)
summary(COL.errW.GM, Hausman=TRUE)
aa <- GMargminImage(COL.errW.GM)
levs <- quantile(aa$z, seq(0, 1, 1/12))
image(aa, breaks=levs, xlab="lambda", ylab="s2")
points(COL.errW.GM$lambda, COL.errW.GM$s2, pch=3, lwd=2)
contour(aa, levels=signif(levs, 4), add=TRUE)
COL.errW.GM1 <- GMerrorsar(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="W"))
summary(COL.errW.GM1)
example(NY_data)
esar1f <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY, family="SAR", method="full")
summary(esar1f)
esar1gm <- GMerrorsar(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
 data=nydata, listw=listw_NY)
summary(esar1gm)
esar1gm1 <- GMerrorsar(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
 data=nydata, listw=listw_NY, method="Nelder-Mead")
summary(esar1gm1)



cleanEx()
nameEx("LR.sarlm")
### * LR.sarlm

flush(stderr()); flush(stdout())

### Name: LR.sarlm
### Title: Likelihood ratio test
### Aliases: LR.sarlm LR1.sarlm Wald1.sarlm Hausman.test Hausman.test.sarlm
###   Hausman.test.gmsar logLik.sarlm
### Keywords: spatial

### ** Examples

example(columbus)
mixed <- lagsarlm(CRIME ~ HOVAL + INC, data=columbus, nb2listw(col.gal.nb),
  type="mixed")
error <- errorsarlm(CRIME ~ HOVAL + INC, data=columbus, nb2listw(col.gal.nb))
LR.sarlm(mixed, error)
Hausman.test(error)



cleanEx()
nameEx("ME")
### * ME

flush(stderr()); flush(stdout())

### Name: ME
### Title: Moran eigenvector GLM filtering
### Aliases: ME print.ME_res fitted.ME_res
### Keywords: spatial

### ** Examples

## Not run: 
##D example(columbus)
##D lmbase <- lm(CRIME ~ INC + HOVAL, data=columbus)
##D lagcol <- SpatialFiltering(CRIME ~ 1, ~ INC + HOVAL, data=columbus,
##D  nb=col.gal.nb, style="W", alpha=0.1, verbose=TRUE)
##D lagcol
##D lmlag <- lm(CRIME ~ INC + HOVAL + fitted(lagcol), data=columbus)
##D anova(lmlag)
##D anova(lmbase, lmlag)
##D set.seed(123)
##D lagcol1 <- ME(CRIME ~ INC + HOVAL, data=columbus, family="gaussian",
##D  listw=nb2listw(col.gal.nb), alpha=0.1, verbose=TRUE)
##D lagcol1
##D lmlag1 <- lm(CRIME ~ INC + HOVAL + fitted(lagcol1), data=columbus)
##D anova(lmlag1)
##D anova(lmbase, lmlag1)
##D set.seed(123)
##D lagcol2 <- ME(CRIME ~ INC + HOVAL, data=columbus, family="gaussian",
##D  listw=nb2listw(col.gal.nb), alpha=0.1, stdev=TRUE, verbose=TRUE)
##D lagcol2
##D lmlag2 <- lm(CRIME ~ INC + HOVAL + fitted(lagcol2), data=columbus)
##D anova(lmlag2)
##D anova(lmbase, lmlag2)
##D example(nc.sids)
##D glmbase <- glm(SID74 ~ 1, data=nc.sids, offset=log(BIR74),
##D  family="poisson")
##D set.seed(123)
##D MEpois1 <- ME(SID74 ~ 1, data=nc.sids, offset=log(BIR74),
##D  family="poisson", listw=nb2listw(ncCR85_nb), alpha=0.2, verbose=TRUE)
##D MEpois1
##D glmME <- glm(SID74 ~ 1 + fitted(MEpois1), data=nc.sids, offset=log(BIR74),
##D  family="poisson")
##D anova(glmME, test="Chisq")
##D anova(glmbase, glmME, test="Chisq")
##D data(hopkins)
##D hopkins_part <- hopkins[21:36,36:21]
##D hopkins_part[which(hopkins_part > 0, arr.ind=TRUE)] <- 1
##D hopkins.rook.nb <- cell2nb(16, 16, type="rook")
##D glmbase <- glm(c(hopkins_part) ~ 1, family="binomial")
##D set.seed(123)
##D MEbinom1 <- ME(c(hopkins_part) ~ 1, family="binomial",
##D  listw=nb2listw(hopkins.rook.nb), alpha=0.2, verbose=TRUE)
##D glmME <- glm(c(hopkins_part) ~ 1 + fitted(MEbinom1), family="binomial")
##D anova(glmME, test="Chisq")
##D anova(glmbase, glmME, test="Chisq")
## End(Not run)



cleanEx()
nameEx("NY_data")
### * NY_data

flush(stderr()); flush(stdout())

### Name: NY_data
### Title: New York leukemia data
### Aliases: NY_data nydata listw_NY
### Keywords: datasets

### ** Examples

## NY leukemia
nydata <- read.dbf(system.file("etc/misc/nydata.dbf", package="spdep")[1])
coordinates(nydata) <- c("X", "Y")
nyadjmat <- as.matrix(read.dbf(system.file("etc/misc/nyadjwts.dbf",
 package="spdep")[1])[-1])
ID <- as.character(names(read.dbf(system.file("etc/misc/nyadjwts.dbf",
 package="spdep")[1]))[-1])
identical(substring(ID, 2, 10), substring(as.character(nydata$AREAKEY), 2, 10))
nyadjlw <- mat2listw(nyadjmat, as.character(nydata$AREAKEY))
listw_NY <- nb2listw(nyadjlw$neighbours, style="B")



cleanEx()
nameEx("SpatialFiltering")
### * SpatialFiltering

flush(stderr()); flush(stdout())

### Name: SpatialFiltering
### Title: Semi-parametric spatial filtering
### Aliases: SpatialFiltering print.SFResult fitted.SFResult
### Keywords: spatial

### ** Examples

example(columbus)
lmbase <- lm(CRIME ~ INC + HOVAL, data=columbus)
sarcol <- SpatialFiltering(CRIME ~ INC + HOVAL, data=columbus,
 nb=col.gal.nb, style="W", ExactEV=TRUE)
sarcol
lmsar <- lm(CRIME ~ INC + HOVAL + fitted(sarcol), data=columbus)
lmsar
anova(lmbase, lmsar)
lm.morantest(lmsar, nb2listw(col.gal.nb))
lagcol <- SpatialFiltering(CRIME ~ 1, ~ INC + HOVAL - 1, data=columbus,
 nb=col.gal.nb, style="W")
lagcol
lmlag <- lm(CRIME ~ INC + HOVAL + fitted(lagcol), data=columbus)
lmlag
anova(lmbase, lmlag)
lm.morantest(lmlag, nb2listw(col.gal.nb))



cleanEx()
nameEx("afcon")
### * afcon

flush(stderr()); flush(stdout())

### Name: afcon
### Title: Spatial patterns of conflict in Africa 1966-78
### Aliases: afcon africa.rook.nb afxy paper.nb
### Keywords: datasets

### ** Examples

data(afcon)
plot(africa.rook.nb, afxy)
plot(diffnb(paper.nb, africa.rook.nb), afxy, col="red", add=TRUE)
text(afxy, labels=attr(africa.rook.nb, "region.id"), pos=4, offset=0.4)
moran.test(afcon$totcon, nb2listw(africa.rook.nb))
moran.test(afcon$totcon, nb2listw(paper.nb))
geary.test(afcon$totcon, nb2listw(paper.nb))



cleanEx()
nameEx("aggregate.nb")
### * aggregate.nb

flush(stderr()); flush(stdout())

### Name: aggregate.nb
### Title: Aggregate a spatial neighbours object
### Aliases: aggregate.nb
### Keywords: spatial

### ** Examples

data(used.cars)
data(state)
cont_st <- match(attr(usa48.nb, "region.id"), state.abb)
cents <- as.matrix(as.data.frame(state.center))[cont_st,]
opar <- par(mfrow=c(2,1))
plot(usa48.nb, cents, xlim=c(-125, -65), ylim=c(25, 50))
IDs <- as.character(state.division[cont_st])
agg_cents <- aggregate(cents, list(IDs), mean)
agg_nb <- aggregate(usa48.nb, IDs)
plot(agg_nb, agg_cents[, 2:3], xlim=c(-125, -65), ylim=c(25, 50))
text(agg_cents[, 2:3], agg_cents[, 1], cex=0.6)
par(opar)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("anova.sarlm")
### * anova.sarlm

flush(stderr()); flush(stdout())

### Name: anova.sarlm
### Title: Comparison of simultaneous autoregressive models
### Aliases: anova.sarlm
### Keywords: spatial

### ** Examples

example(columbus)
lm.mod <- lm(CRIME ~ HOVAL + INC, data=columbus)
lag <- lagsarlm(CRIME ~ HOVAL + INC, data=columbus, nb2listw(col.gal.nb))
mixed <- lagsarlm(CRIME ~ HOVAL + INC, data=columbus, nb2listw(col.gal.nb),
  type="mixed")
error <- errorsarlm(CRIME ~ HOVAL + INC, data=columbus, nb2listw(col.gal.nb))
LR.sarlm(mixed, error)
anova(lag, lm.mod)
anova(lag, error, mixed)
AIC(lag, error, mixed)



cleanEx()
nameEx("aple")
### * aple

flush(stderr()); flush(stdout())

### Name: aple
### Title: Approximate profile-likelihood estimator (APLE)
### Aliases: aple
### Keywords: spatial

### ** Examples

example(wheat)
nbr1 <- poly2nb(wheat, queen=FALSE)
nbrl <- nblag(nbr1, 2)
nbr12 <- nblag_cumul(nbrl)
cms0 <- with(as(wheat, "data.frame"), tapply(yield, c, median))
cms1 <- c(model.matrix(~ factor(c) -1, data=wheat) %*% cms0)
wheat$yield_detrend <- wheat$yield - cms1
isTRUE(all.equal(c(with(as(wheat, "data.frame"),
 tapply(yield_detrend, c, median))), rep(0.0, 25),
 check.attributes=FALSE))
moran.test(wheat$yield_detrend, nb2listw(nbr12, style="W"))
aple(as.vector(scale(wheat$yield_detrend, scale=FALSE)), nb2listw(nbr12, style="W"))
errorsarlm(yield_detrend ~ 1, wheat, nb2listw(nbr12, style="W"))



cleanEx()
nameEx("aple.mc")
### * aple.mc

flush(stderr()); flush(stdout())

### Name: aple.mc
### Title: Approximate profile-likelihood estimator (APLE) permutation test
### Aliases: aple.mc
### Keywords: spatial

### ** Examples

## Not run: 
##D example(aple)
##D boot_out <- aple.mc(as.vector(scale(wheat$yield_detrend, scale=FALSE)),
##D  nb2listw(nbr12, style="W"), nsim=500)
##D plot(boot_out)
##D boot_out
## End(Not run)



cleanEx()
nameEx("aple.plot")
### * aple.plot

flush(stderr()); flush(stdout())

### Name: aple.plot
### Title: Approximate profile-likelihood estimator (APLE) scatterplot
### Aliases: aple.plot localAple
### Keywords: spatial

### ** Examples

## Not run: 
##D example(aple)
##D plt_out <- aple.plot(as.vector(scale(wheat$yield_detrend, scale=FALSE)),
##D  nb2listw(nbr12, style="W"), cex=0.6)
##D crossprod(plt_out$Y, plt_out$X)/crossprod(plt_out$X)
##D lm_obj <- lm(Y ~ X, plt_out)
##D abline(lm_obj)
##D abline(v=0, h=0, lty=2)
##D zz <- summary(influence.measures(lm_obj))
##D infl <- as.integer(rownames(zz))
##D points(plt_out$X[infl], plt_out$Y[infl], pch=3, cex=0.6, col="red")
##D wheat$localAple <- localAple(as.vector(scale(wheat$yield_detrend, scale=FALSE)),
##D  nb2listw(nbr12, style="W"))
##D mean(wheat$localAple)
##D hist(wheat$localAple)
##D spl <- list("sp.text", coordinates(wheat)[infl,], rep("*", length(infl)))
##D spplot(wheat, "localAple", sp.layout=spl)
## End(Not run)



cleanEx()
nameEx("as_dgRMatrix_listw")
### * as_dgRMatrix_listw

flush(stderr()); flush(stdout())

### Name: as_dgRMatrix_listw
### Title: Interface between Matrix class objects and weights lists
### Aliases: as_dgRMatrix_listw as_dsTMatrix_listw as_dsCMatrix_I
###   as_dsCMatrix_IrW Jacobian_W
### Keywords: spatial

### ** Examples

example(NY_data)
W <- as_dsTMatrix_listw(listw_NY)
I <- as_dsCMatrix_I(dim(W)[1])
W <- as(W, "CsparseMatrix")
rho <- 0.1
c(determinant(I - rho * W, logarithm=TRUE)$modulus)
sum(log(1 - rho * eigenw(listw_NY)))
n <- dim(W)[1]
nW <- - W
nChol <- Cholesky(nW, Imult=8)
.f <- if(package_version(packageDescription("Matrix")$Version) >
           "0.999375-30") 2 else 1
n * log(rho) + (.f * c(determinant(update(nChol, nW, 1/rho))$modulus))
rho <- seq(0.01, 0.1, 0.01)
n * log(rho) + Matrix:::ldetL2up(nChol, nW, 1/rho)



cleanEx()
nameEx("auckland")
### * auckland

flush(stderr()); flush(stdout())

### Name: auckland
### Title: Marshall's infant mortality in Auckland dataset
### Aliases: auckland auckland.nb auckpolys
### Keywords: datasets

### ** Examples

auckland <- readShapePoly(system.file("etc/shapes/auckland.shp",
 package="spdep")[1])
auckland.nb <- poly2nb(auckland)



cleanEx()
nameEx("autocov_dist")
### * autocov_dist

flush(stderr()); flush(stdout())

### Name: autocov_dist
### Title: Distance-weighted autocovariate
### Aliases: autocov_dist
### Keywords: spatial

### ** Examples

example(columbus)
xy <- cbind(columbus$X, columbus$Y)
ac1a <- autocov_dist(columbus$CRIME, xy, nbs=10, style="W",
 type="one")
acinva <- autocov_dist(columbus$CRIME, xy, nbs=10, style="W",
 type="inverse")
acinv2a <- autocov_dist(columbus$CRIME, xy, nbs=10, style="W",
 type="inverse.squared")

plot(ac1a ~ columbus$CRIME, pch=16, asp=1)
points(acinva ~ columbus$CRIME, pch=16, col="red")
points(acinv2a ~ columbus$CRIME, pch=16, col="blue")
abline(0,1)

nb <- dnearneigh(xy, 0, 10)
lw <- nb2listw(nb, style="W")
ac1b <- lag(lw, columbus$CRIME)
all.equal(ac1b, ac1a)

nbd <- nbdists(nb, xy)
gl <- lapply(nbd, function(x) 1/x)
lw <- nb2listw(nb, glist=gl)
acinvb <- lag(lw, columbus$CRIME)
all.equal(acinvb, acinva)

gl2 <- lapply(nbd, function(x) 1/(x^2))
lw <- nb2listw(nb, glist=gl2)
acinv2b <- lag(lw, columbus$CRIME)
all.equal(acinv2b, acinv2a)

glm(CRIME ~ HOVAL + ac1b, data=columbus, family="gaussian")
spautolm(columbus$CRIME ~ HOVAL, data=columbus,
 listw=nb2listw(nb, style="W"))

xy <- SpatialPoints(xy)
acinva <- autocov_dist(columbus$CRIME, xy, nbs=10, style="W",
 type="inverse")
nb <- dnearneigh(xy, 0, 10)
nbd <- nbdists(nb, xy)
gl <- lapply(nbd, function(x) 1/x)
lw <- nb2listw(nb, glist=gl)
acinvb <- lag(lw, columbus$CRIME)
all.equal(acinvb, acinva)




cleanEx()
nameEx("baltimore")
### * baltimore

flush(stderr()); flush(stdout())

### Name: baltimore
### Title: House sales prices, Baltimore, MD 1978
### Aliases: baltimore
### Keywords: datasets

### ** Examples

data(baltimore)
## maybe str(baltimore) ; plot(baltimore) ...



cleanEx()
nameEx("bhicv")
### * bhicv

flush(stderr()); flush(stdout())

### Name: bhicv
### Title: Data set with 4 life condition indices of Belo Horizonte region
### Aliases: bhicv
### Keywords: data

### ** Examples

### see example in 'skater' function help



cleanEx()
nameEx("boston")
### * boston

flush(stderr()); flush(stdout())

### Name: boston
### Title: Corrected Boston Housing Data
### Aliases: boston.c boston boston.soi boston.utm
### Keywords: datasets

### ** Examples

data(boston)
hr0 <- lm(log(MEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + I(RM^2) +
 AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT), data=boston.c)
summary(hr0)
logLik(hr0)
gp0 <- lm(log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + I(RM^2) +
 AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT), data=boston.c)
summary(gp0)
logLik(gp0)
lm.morantest(hr0, nb2listw(boston.soi))
## Not run: 
##D gp1 <- errorsarlm(log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2)
##D  + I(RM^2) +  AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT),
##D  data=boston.c, nb2listw(boston.soi), method="Matrix", 
##D  control=list(tol.opt = .Machine$double.eps^(1/4)))
##D summary(gp1)
##D gp2 <- lagsarlm(log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + I(RM^2)
##D  +  AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT),
##D  data=boston.c, nb2listw(boston.soi), method="Matrix")
##D summary(gp2)
## End(Not run)




cleanEx()
nameEx("bptest.sarlm")
### * bptest.sarlm

flush(stderr()); flush(stdout())

### Name: bptest.sarlm
### Title: Breusch-Pagan test for spatial models
### Aliases: bptest.sarlm
### Keywords: spatial

### ** Examples

example(columbus)
error.col <- errorsarlm(CRIME ~ HOVAL + INC, data=columbus,
 nb2listw(col.gal.nb))
bptest.sarlm(error.col)
bptest.sarlm(error.col, studentize=FALSE)
## Not run: 
##D lm.target <- lm(error.col$tary ~ error.col$tarX - 1)
##D if (require(lmtest) && require(sandwich)) {
##D   coeftest(lm.target, vcov=vcovHC(lm.target, type="HC0"), df=Inf)
##D }
## End(Not run)



cleanEx()
nameEx("card")
### * card

flush(stderr()); flush(stdout())

### Name: card
### Title: Cardinalities for neighbours lists
### Aliases: card
### Keywords: spatial

### ** Examples

example(columbus)
table(card(col.gal.nb))



cleanEx()
nameEx("cell2nb")
### * cell2nb

flush(stderr()); flush(stdout())

### Name: cell2nb
### Title: Generate neighbours list for grid cells
### Aliases: cell2nb mrc2vi rookcell queencell vi2mrc
### Keywords: spatial

### ** Examples

nb7rt <- cell2nb(7, 7)
summary(nb7rt)
xyc <- attr(nb7rt, "region.id")
xy <- matrix(as.integer(unlist(strsplit(xyc, ":"))), ncol=2, byrow=TRUE)
plot(nb7rt, xy)
nb7rt <- cell2nb(7, 7, torus=TRUE)
summary(nb7rt)



cleanEx()
nameEx("choynowski")
### * choynowski

flush(stderr()); flush(stdout())

### Name: choynowski
### Title: Choynowski probability map values
### Aliases: choynowski
### Keywords: spatial

### ** Examples

example(auckland)
res <- choynowski(auckland$M77_85, 9*auckland$Und5_81)
resl <- choynowski(auckland$M77_85, 9*auckland$Und5_81, legacy=TRUE)
all.equal(res, resl)
rt <- sum(auckland$M77_85)/sum(9*auckland$Und5_81)
ch_ppois_pmap <- numeric(length(auckland$Und5_81))
side <- c("greater", "less")
for (i in seq(along=ch_ppois_pmap)) {
  ch_ppois_pmap[i] <- poisson.test(auckland$M77_85[i], r=rt,
    T=(9*auckland$Und5_81[i]), alternative=side[(res$type[i]+1)])$p.value
}
all.equal(ch_ppois_pmap, res$pmap)

res1 <- probmap(auckland$M77_85, 9*auckland$Und5_81)
table(abs(res$pmap - res1$pmap) < 0.00001, res$type)
lt005 <- (res$pmap < 0.05) & (res$type)
ge005 <- (res$pmap < 0.05) & (!res$type)
cols <- rep("white", length(lt005))
cols[lt005] <- grey(2/7)
cols[ge005] <- grey(5/7)
plot(auckland, col=cols) 
legend("bottomleft", fill=grey(c(2,5)/7), legend=c("low", "high"), bty="n")



cleanEx()
nameEx("columbus")
### * columbus

flush(stderr()); flush(stdout())

### Name: columbus
### Title: Columbus OH spatial analysis data set
### Aliases: columbus col.gal.nb coords polys bbs
### Keywords: datasets

### ** Examples

columbus <- readShapePoly(system.file("etc/shapes/columbus.shp",
 package="spdep")[1])
col.gal.nb <- read.gal(system.file("etc/weights/columbus.gal",
 package="spdep")[1])



cleanEx()
nameEx("compon")
### * compon

flush(stderr()); flush(stdout())

### Name: Graph Components
### Title: Depth First Search on Neighbor Lists
### Aliases: n.comp.nb
### Keywords: spatial

### ** Examples

example(columbus)
coords <- coordinates(columbus)
plot(col.gal.nb, coords, col="grey")
col2 <- droplinks(col.gal.nb, 21)
plot(col2, coords, add=TRUE)
res <- n.comp.nb(col2)
table(res$comp.id)
points(coords, col=res$comp.id, pch=16)



cleanEx()
nameEx("diffnb")
### * diffnb

flush(stderr()); flush(stdout())

### Name: diffnb
### Title: Differences between neighbours lists
### Aliases: diffnb
### Keywords: spatial

### ** Examples

example(columbus)
coords <- coordinates(columbus)
rn <- sapply(slot(columbus, "polygons"), function(x) slot(x, "ID"))
knn1 <- knearneigh(coords, 1)
knn2 <- knearneigh(coords, 2)
nb1 <- knn2nb(knn1, row.names=rn)
nb2 <- knn2nb(knn2, row.names=rn)
diffs <- diffnb(nb2, nb1)
plot(columbus, border="grey")
plot(nb1, coords, add=TRUE)
plot(diffs, coords, add=TRUE, col="red", lty=2)
title(main="Plot of first (black) and second (red)\nnearest neighbours")



cleanEx()
nameEx("dnearneigh")
### * dnearneigh

flush(stderr()); flush(stdout())

### Name: dnearneigh
### Title: Neighbourhood contiguity by distance
### Aliases: dnearneigh
### Keywords: spatial

### ** Examples

example(columbus)
coords <- coordinates(columbus)
rn <- sapply(slot(columbus, "polygons"), function(x) slot(x, "ID"))
k1 <- knn2nb(knearneigh(coords))
all.linked <- max(unlist(nbdists(k1, coords)))
col.nb.0.all <- dnearneigh(coords, 0, all.linked, row.names=rn)
summary(col.nb.0.all, coords)
plot(columbus, border="grey")
plot(col.nb.0.all, coords, add=TRUE)
title(main=paste("Distance based neighbours 0-",  format(all.linked),
 " distance units", sep=""))
data(state)
us48.fipsno <- read.geoda(system.file("etc/weights/us48.txt",
 package="spdep")[1])
if (as.numeric(paste(version$major, version$minor, sep="")) < 19) {
 m50.48 <- match(us48.fipsno$"State.name", state.name)
} else {
 m50.48 <- match(us48.fipsno$"State_name", state.name)
}
xy <- as.matrix(as.data.frame(state.center))[m50.48,]
llk1 <- knn2nb(knearneigh(xy, k=1, longlat=FALSE))
all.linked <- max(unlist(nbdists(llk1, xy, longlat=FALSE)))
ll.nb <- dnearneigh(xy, 0, all.linked, longlat=FALSE)
summary(ll.nb, xy, longlat=TRUE, scale=0.5)
gck1 <- knn2nb(knearneigh(xy, k=1, longlat=TRUE))
all.linked <- max(unlist(nbdists(gck1, xy, longlat=TRUE)))
gc.nb <- dnearneigh(xy, 0, all.linked, longlat=TRUE)
summary(gc.nb, xy, longlat=TRUE, scale=0.5)
plot(ll.nb, xy)
plot(diffnb(ll.nb, gc.nb), xy, add=TRUE, col="red", lty=2)
title(main="Differences between Euclidean and Great Circle neighbours")

xy1 <- SpatialPoints((as.data.frame(state.center))[m50.48,],
  proj4string=CRS("+proj=longlat"))
gck1a <- knn2nb(knearneigh(xy1, k=1))
all.linked <- max(unlist(nbdists(gck1a, xy1)))
gc.nb <- dnearneigh(xy1, 0, all.linked)
summary(gc.nb, xy1, scale=0.5)



cleanEx()
nameEx("do_ldet")
### * do_ldet

flush(stderr()); flush(stdout())

### Name: do_ldet
### Title: Spatial regression model Jacobian computations
### Aliases: do_ldet eigen_setup mcdet_setup cheb_setup spam_setup
###   spam_update_setup Matrix_setup Matrix_J_setup LU_setup moments_setup
### Keywords: spatial

### ** Examples

data(boston)
lw <- nb2listw(boston.soi)
can.sim <- spdep:::can.be.simmed(lw)
env <- new.env(parent=globalenv())
assign("listw", lw, envir=env)
assign("can.sim", can.sim, envir=env)
assign("similar", FALSE, envir=env)
assign("verbose", FALSE, envir=env)
assign("family", "SAR", envir=env)
eigen_setup(env)
get("similar", envir=env)
do_ldet(0.5, env)
rm(env)
env <- new.env(parent=globalenv())
assign("listw", lw, envir=env)
assign("can.sim", can.sim, envir=env)
assign("similar", FALSE, envir=env)
assign("family", "SAR", envir=env)
assign("n", length(boston.soi), envir=env)
Matrix_setup(env, Imult=2, super=FALSE)
get("similar", envir=env)
do_ldet(0.5, env)
rm(env)
env <- new.env(parent=globalenv())
assign("listw", lw, envir=env)
assign("n", length(boston.soi), envir=env)
assign("can.sim", can.sim, envir=env)
assign("similar", FALSE, envir=env)
assign("family", "SAR", envir=env)
spam_setup(env)
get("similar", envir=env)
do_ldet(0.5, env)
rm(env)
env <- new.env(parent=globalenv())
assign("listw", lw, envir=env)
assign("n", length(boston.soi), envir=env)
assign("similar", FALSE, envir=env)
assign("family", "SAR", envir=env)
LU_setup(env)
get("similar", envir=env)
do_ldet(0.5, env)
rm(env)
env <- new.env(parent=globalenv())
assign("listw", lw, envir=env)
assign("similar", FALSE, envir=env)
assign("family", "SAR", envir=env)
cheb_setup(env, q=5)
get("similar", envir=env)
do_ldet(0.5, env)
rm(env)
env <- new.env(parent=globalenv())
assign("listw", lw, envir=env)
assign("n", length(boston.soi), envir=env)
assign("similar", FALSE, envir=env)
assign("family", "SAR", envir=env)
set.seed(12345)
mcdet_setup(env, p=16, m=30)
get("similar", envir=env)
do_ldet(0.5, env)
rm(env)



cleanEx()
nameEx("droplinks")
### * droplinks

flush(stderr()); flush(stdout())

### Name: droplinks
### Title: Drop links in a neighbours list
### Aliases: droplinks
### Keywords: spatial

### ** Examples

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
    res <- diag(solve(t(mat) %*% mat))
    mns[i,j] <- mean(res)
    mxs[i,j] <- max(res)
  }
}
print(mns)
print(mxs)





cleanEx()
nameEx("edit.nb")
### * edit.nb

flush(stderr()); flush(stdout())

### Name: edit.nb
### Title: Interactive editing of neighbours lists
### Aliases: edit.nb
### Keywords: spatial

### ** Examples

## Not run: 
##D data(columbus)
##D class(polys)
##D nnb <- edit.nb(col.gal.nb, coords, polys)
##D example(columbus)
##D class(columbus)
##D nnb1 <- edit.nb(col.gal.nb, polys=columbus)
## End(Not run)


cleanEx()
nameEx("eigenw")
### * eigenw

flush(stderr()); flush(stdout())

### Name: eigenw
### Title: Spatial weights matrix eigenvalues
### Aliases: eigenw
### Keywords: spatial

### ** Examples

data(oldcol)
W.eig <- eigenw(nb2listw(COL.nb, style="W"))
1/range(W.eig)
S.eig <- eigenw(nb2listw(COL.nb, style="S"))
1/range(S.eig)
B.eig <- eigenw(nb2listw(COL.nb, style="B"))
1/range(B.eig)
# cases for intrinsically asymmetric weights
crds <- cbind(COL.OLD$X, COL.OLD$Y)
k6 <- knn2nb(knearneigh(crds, k=6))
is.symmetric.nb(k6)
k6eig <- eigenw(nb2listw(k6, style="W"))
is.complex(k6eig)
rho <- 0.5
Jc <- sum(log(1 - rho * k6eig))
# complex eigenvalue Jacobian
Jc
W <- as(as_dgRMatrix_listw(nb2listw(k6, style="W")), "CsparseMatrix")
I <- diag(length(k6))
Jl <- sum(log(abs(diag(slot(lu(I - rho * W), "U")))))
# LU Jacobian equals complex eigenvalue Jacobian
Jl
all.equal(Re(Jc), Jl)
# wrong value if only real part used
Jr <- sum(log(1 - rho * Re(k6eig)))
Jr
all.equal(Jr, Jl)
# construction of Jacobian from complex conjugate pairs (Jan Hauke)
Rev <- Re(k6eig)[which(Im(k6eig) == 0)]
# real eigenvalues
Cev <- k6eig[which(Im(k6eig) != 0)]
pCev <- Cev[Im(Cev) > 0]
# separate complex conjugate pairs
RpCev <- Re(pCev)
IpCev <- Im(pCev)
# reassemble Jacobian
Jc1 <- sum(log(1 - rho*Rev)) + sum(log((1 - rho * RpCev)^2 + (rho^2)*(IpCev^2)))
all.equal(Re(Jc), Jc1)
# impact of omitted complex part term in real part only Jacobian
Jc2 <- sum(log(1 - rho*Rev)) + sum(log((1 - rho * RpCev)^2))
all.equal(Jr, Jc2)
# trace of asymmetric (WW) and crossprod of complex eigenvalues for APLE
sum(diag(W %*% W))
crossprod(k6eig)



cleanEx()
nameEx("eire")
### * eire

flush(stderr()); flush(stdout())

### Name: eire
### Title: Eire data sets
### Aliases: eire eire.df eire.polys.utm eire.coords.utm eire.nb
### Keywords: datasets

### ** Examples

eire <- readShapePoly(system.file("etc/shapes/eire.shp", package="spdep")[1],
  ID="names", proj4string=CRS("+proj=utm +zone=30 +units=km"))
eire.nb <- poly2nb(eire)
#data(eire)
summary(eire$A)
brks <- round(fivenum(eire$A), digits=2)
cols <- rev(heat.colors(4))
plot(eire, col=cols[findInterval(eire$A, brks, all.inside=TRUE)])
title(main="Percentage with blood group A in Eire")
legend(x=c(-50, 70), y=c(6120, 6050), leglabs(brks), fill=cols, bty="n")
plot(eire)
plot(eire.nb, coordinates(eire), add=TRUE)
lA <- lag.listw(nb2listw(eire.nb), eire$A)
summary(lA)
moran.test(eire$A, nb2listw(eire.nb))
geary.test(eire$A, nb2listw(eire.nb))
cor(lA, eire$A)
moran.plot(eire$A, nb2listw(eire.nb),
 labels=eire$names)
A.lm <- lm(A ~ towns + pale, data=eire)
summary(A.lm)
res <- residuals(A.lm)
brks <- c(min(res),-2,-1,0,1,2,max(res))
cols <- rev(cm.colors(6))
plot(eire, col=cols[findInterval(res, brks, all.inside=TRUE)])
title(main="Regression residuals")
legend(x=c(-50, 70), y=c(6120, 6050), legend=leglabs(brks), fill=cols,
  bty="n")
lm.morantest(A.lm, nb2listw(eire.nb))
lm.morantest.sad(A.lm, nb2listw(eire.nb))
lm.LMtests(A.lm, nb2listw(eire.nb), test="LMerr")
brks <- round(fivenum(eire$OWNCONS), digits=2)
cols <- grey(4:1/5)
plot(eire, col=cols[findInterval(eire$OWNCONS, brks, all.inside=TRUE)])
title(main="Percentage own consumption of agricultural produce")
legend(x=c(-50, 70), y=c(6120, 6050), legend=leglabs(brks),
  fill=cols, bty="n")
moran.plot(eire$OWNCONS, nb2listw(eire.nb))
moran.test(eire$OWNCONS, nb2listw(eire.nb))
e.lm <- lm(OWNCONS ~ ROADACC, data=eire)
res <- residuals(e.lm)
brks <- c(min(res),-2,-1,0,1,2,max(res))
cols <- rev(cm.colors(6))
plot(eire, col=cols[findInterval(res, brks, all.inside=TRUE)])
title(main="Regression residuals")
legend(x=c(-50, 70), y=c(6120, 6050), legend=leglabs(brks), fill=cm.colors(6),
  bty="n")
lm.morantest(e.lm, nb2listw(eire.nb))
lm.morantest.sad(e.lm, nb2listw(eire.nb))
lm.LMtests(e.lm, nb2listw(eire.nb), test="LMerr")
print(localmoran.sad(e.lm, eire.nb, select=1:length(slot(eire, "polygons"))))



cleanEx()
nameEx("elect80")
### * elect80

flush(stderr()); flush(stdout())

### Name: elect80
### Title: 1980 Presidential election results
### Aliases: elect80 elect80_lw k4 dll e80_queen
### Keywords: datasets

### ** Examples

data(elect80)



cleanEx()
nameEx("errorsarlm")
### * errorsarlm

flush(stderr()); flush(stdout())

### Name: errorsarlm
### Title: Spatial simultaneous autoregressive error model estimation
### Aliases: errorsarlm
### Keywords: spatial

### ** Examples

data(oldcol)
COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="W"), method="eigen", quiet=FALSE)
summary(COL.errW.eig, correlation=TRUE)
COL.errB.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="B"), method="eigen", quiet=FALSE)
summary(COL.errB.eig, correlation=TRUE)
W <- as(as_dgRMatrix_listw(nb2listw(COL.nb)), "CsparseMatrix")
trMatc <- trW(W, type="mult")
COL.errW.M <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="W"), method="Matrix", quiet=FALSE, trs=trMatc)
summary(COL.errW.M, correlation=TRUE)
NA.COL.OLD <- COL.OLD
NA.COL.OLD$CRIME[20:25] <- NA
COL.err.NA <- errorsarlm(CRIME ~ INC + HOVAL, data=NA.COL.OLD,
 nb2listw(COL.nb), na.action=na.exclude)
COL.err.NA$na.action
COL.err.NA
resid(COL.err.NA)
system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="W"), method="eigen"))
system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="W"), method="eigen", control=list(LAPACK=TRUE)))
system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="W"), method="eigen", control=list(compiled_sse=TRUE)))
system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="W"), method="Matrix_J", control=list(super=TRUE)))
system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="W"), method="Matrix_J", control=list(super=FALSE)))
system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="W"), method="Matrix_J", control=list(super=as.logical(NA))))
system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="W"), method="Matrix", control=list(super=TRUE)))
system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="W"), method="Matrix", control=list(super=FALSE)))
system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="W"), method="Matrix", control=list(super=as.logical(NA))))
system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="W"), method="spam", control=list(spamPivot="MMD")))
system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="W"), method="spam", control=list(spamPivot="RCM")))
system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="W"), method="spam_update", control=list(spamPivot="MMD")))
system.time(COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="W"), method="spam_update", control=list(spamPivot="RCM")))
COL.merrW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="W"), method="eigen", etype="emixed")
summary(COL.merrW.eig, correlation=TRUE)



cleanEx()
nameEx("geary")
### * geary

flush(stderr()); flush(stdout())

### Name: geary
### Title: Compute Geary's C
### Aliases: geary
### Keywords: spatial

### ** Examples

data(oldcol)
col.W <- nb2listw(COL.nb, style="W")
str(geary(COL.OLD$CRIME, col.W, length(COL.nb), length(COL.nb)-1,
 Szero(col.W)))



cleanEx()
nameEx("geary.mc")
### * geary.mc

flush(stderr()); flush(stdout())

### Name: geary.mc
### Title: Permutation test for Geary's C statistic
### Aliases: geary.mc
### Keywords: spatial

### ** Examples

data(oldcol)
sim1 <- geary.mc(COL.OLD$CRIME, nb2listw(COL.nb, style="W"),
 nsim=99, alternative="less")
sim1
mean(sim1$res)
var(sim1$res)
summary(sim1$res)
colold.lags <- nblag(COL.nb, 3)
sim2 <- geary.mc(COL.OLD$CRIME, nb2listw(colold.lags[[2]],
 style="W"), nsim=99)
sim2
summary(sim2$res)
sim3 <- geary.mc(COL.OLD$CRIME, nb2listw(colold.lags[[3]],
 style="W"), nsim=99)
sim3
summary(sim3$res)



cleanEx()
nameEx("geary.test")
### * geary.test

flush(stderr()); flush(stdout())

### Name: geary.test
### Title: Geary's C test for spatial autocorrelation
### Aliases: geary.test
### Keywords: spatial

### ** Examples

data(oldcol)
geary.test(COL.OLD$CRIME, nb2listw(COL.nb, style="W"))
geary.test(COL.OLD$CRIME, nb2listw(COL.nb, style="W"),
 randomisation=FALSE)
colold.lags <- nblag(COL.nb, 3)
geary.test(COL.OLD$CRIME, nb2listw(colold.lags[[2]],
 style="W"))
geary.test(COL.OLD$CRIME, nb2listw(colold.lags[[3]],
 style="W"), alternative="greater")
print(is.symmetric.nb(COL.nb))
coords.OLD <- cbind(COL.OLD$X, COL.OLD$Y)
COL.k4.nb <- knn2nb(knearneigh(coords.OLD, 4))
print(is.symmetric.nb(COL.k4.nb))
geary.test(COL.OLD$CRIME, nb2listw(COL.k4.nb, style="W"))
geary.test(COL.OLD$CRIME, nb2listw(COL.k4.nb, style="W"),
 randomisation=FALSE)
cat("Note non-symmetric weights matrix - use listw2U()\n")
geary.test(COL.OLD$CRIME, listw2U(nb2listw(COL.k4.nb,
 style="W")))
geary.test(COL.OLD$CRIME, listw2U(nb2listw(COL.k4.nb,
 style="W")), randomisation=FALSE)



cleanEx()
nameEx("getisord")
### * getisord

flush(stderr()); flush(stdout())

### Name: getisord
### Title: Getis-Ord remote sensing example data
### Aliases: getisord x y xyz
### Keywords: datasets

### ** Examples

data(getisord)
image(x, y, t(matrix(xyz$val, nrow=16, ncol=16, byrow=TRUE)), asp=1)
text(xyz$x, xyz$y, xyz$val, cex=0.7)
polygon(c(195,225,225,195), c(195,195,225,225), lwd=2)
title(main="Getis-Ord 1996 remote sensing data")



cleanEx()
nameEx("globalG.test")
### * globalG.test

flush(stderr()); flush(stdout())

### Name: globalG.test
### Title: Global G test for spatial autocorrelation
### Aliases: globalG.test
### Keywords: spatial

### ** Examples

example(nc.sids)
sidsrate79 <- (1000*nc.sids$SID79)/nc.sids$BIR79
dists <- c(10, 20, 30, 33, 40, 50, 60, 70, 80, 90, 100)
ndists <- length(dists)
ZG <- numeric(length=ndists)
milesxy <- cbind(nc.sids$east, nc.sids$north)
for (i in 1:ndists) {
  thisnb <- dnearneigh(milesxy, 0, dists[i])
  thislw <- nb2listw(thisnb, style="B", zero.policy=TRUE)
  ZG[i] <- globalG.test(sidsrate79, thislw, zero.policy=TRUE)$statistic
}
cbind(dists, ZG)



cleanEx()
nameEx("graphneigh")
### * graphneigh

flush(stderr()); flush(stdout())

### Name: graphneigh
### Title: Graph based spatial weights
### Aliases: gabrielneigh relativeneigh soi.graph plot.Gabriel
###   plot.relative graph2nb
### Keywords: spatial

### ** Examples

example(columbus)
coords <- coordinates(columbus)
par(mfrow=c(2,2))
col.tri.nb<-tri2nb(coords)
col.gab.nb<-graph2nb(gabrielneigh(coords), sym=TRUE)
col.rel.nb<- graph2nb(relativeneigh(coords), sym=TRUE)
col.soi.nb<- graph2nb(soi.graph(col.tri.nb,coords), sym=TRUE)
plot(columbus, border="grey")
plot(col.tri.nb,coords,add=TRUE)
title(main="Delaunay Triangulation")
plot(columbus, border="grey")
plot(col.gab.nb, coords, add=TRUE)
title(main="Gabriel Graph")
plot(columbus, border="grey")
plot(col.rel.nb, coords, add=TRUE)
title(main="Relative Neighbor Graph")
plot(columbus, border="grey")
plot(col.soi.nb, coords, add=TRUE)
title(main="Sphere of Influence Graph")
par(mfrow=c(1,1))
dx <- rep(0.25*0:4,5)
dy <- c(rep(0,5),rep(0.25,5),rep(0.5,5), rep(0.75,5),rep(1,5))
m <- cbind(c(dx, dx, 3+dx, 3+dx), c(dy, 3+dy, dy, 3+dy))
try(res <- gabrielneigh(m))
res <- gabrielneigh(m, nnmult=4)
summary(graph2nb(res))



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("gstsls")
### * gstsls

flush(stderr()); flush(stdout())

### Name: gstsls
### Title: Spatial simultaneous autoregressive SAC model estimation by GMM
### Aliases: gstsls
### Keywords: spatial

### ** Examples

data(oldcol)
COL.errW.GM <- gstsls(CRIME ~ INC + HOVAL, data=COL.OLD, nb2listw(COL.nb, style="W"))
summary(COL.errW.GM)
aa <- GMargminImage(COL.errW.GM)
levs <- quantile(aa$z, seq(0, 1, 1/12))
image(aa, breaks=levs, xlab="lambda", ylab="s2")
points(COL.errW.GM$lambda, COL.errW.GM$s2, pch=3, lwd=2)
contour(aa, levels=signif(levs, 4), add=TRUE)
COL.errW.GM <- gstsls(CRIME ~ INC + HOVAL, data=COL.OLD, nb2listw(COL.nb, style="W"), scaleU=TRUE)
summary(COL.errW.GM)
listw <- nb2listw(COL.nb)
W <- as(as_dgRMatrix_listw(listw), "CsparseMatrix")
trMat <- trW(W, type="mult")
impacts(COL.errW.GM, tr=trMat)



cleanEx()
nameEx("hopkins")
### * hopkins

flush(stderr()); flush(stdout())

### Name: hopkins
### Title: Hopkins burnt savanna herb remains
### Aliases: hopkins
### Keywords: datasets

### ** Examples

data(hopkins)
image(1:32, 1:32, hopkins[5:36,36:5], breaks=c(-0.5, 3.5, 20),
 col=c("white", "black"))
box()



cleanEx()
nameEx("house")
### * house

flush(stderr()); flush(stdout())

### Name: house
### Title: Lucas county OH housing
### Aliases: house LO_nb trMat
### Keywords: datasets

### ** Examples

## Not run: 
##D house <- read.table("house.dat", header=FALSE)
##D names(house) <- c("price", "yrbuilt", "stories", "TLA", "wall", "beds",
##D   "baths", "halfbaths", "frontage", "depth", "garage", "garagesqft", "rooms",
##D   "lotsize", "sdate", "avalue", "long", "lat", "s1993", "s1994", "s1995",
##D   "s1996", "s1997", "s1998")
##D house$syear <- 1992 + house$s1993 + 2*house$s1994 + 3*house$s1995 +
##D 4*house$s1996 + 5*house$s1997 + 6*house$s1998
##D house$syear <- factor(house$syear)
##D house$age <- (1999 - house$yrbuilt)/100
##D house$stories <- factor(house$stories, levels=1:7, labels=c("one",
##D  "bilevel", "multilvl", "one+half", "two", "two+half", "three"))
##D house$wall <- factor(house$wall, levels=1:7, labels=c("stucdrvt",
##D  "ccbtile", "metlvnyl", "brick", "stone", "wood", "partbrk"))
##D house$garage <- factor(house$garage, levels=0:4, labels=c("no garage",
##D  "basement", "attached", "detached", "carport"))
##D library(sp)
##D coordinates(house) <- c("long", "lat")
##D proj4string(house) <- CRS("+proj=longlat")
##D library(rgdal)
##D house <- spTransform(house, CRS("+init=epsg:2834"))
##D library(spdep)
##D LO_nb <- graph2nb(soi.graph(tri2nb(coordinates(house)), coordinates(house)))
##D W <- as(as_dgRMatrix_listw(nb2listw(LO_nb)), "CsparseMatrix")
##D trMat <- trW(W, type="mult")
## End(Not run)
data(house)
## maybe str(house) ; plot(house) ...



cleanEx()
nameEx("huddersfield")
### * huddersfield

flush(stderr()); flush(stdout())

### Name: huddersfield
### Title: Prevalence of respiratory symptoms
### Aliases: huddersfield
### Keywords: datasets

### ** Examples

data(huddersfield)
str(huddersfield)



cleanEx()
nameEx("impacts.sarlm")
### * impacts.sarlm

flush(stderr()); flush(stdout())

### Name: impacts
### Title: Impacts in spatial lag models
### Aliases: impacts impacts.sarlm impacts.stsls impacts.gmsar
###   plot.lagImpact print.lagImpact summary.lagImpact
###   print.summary.lagImpact HPDinterval.lagImpact
### Keywords: spatial

### ** Examples

example(columbus)
listw <- nb2listw(col.gal.nb)
lobj <- lagsarlm(CRIME ~ INC + HOVAL, columbus, listw)
summary(lobj)
mobj <- lagsarlm(CRIME ~ INC + HOVAL, columbus, listw, type="mixed")
summary(mobj)
W <- as(as_dgRMatrix_listw(listw), "CsparseMatrix")
trMatc <- trW(W, type="mult")
trMC <- trW(W, type="MC")
impacts(lobj, listw=listw)
impacts(lobj, tr=trMatc)
impacts(lobj, tr=trMC)
lobj1 <- stsls(CRIME ~ INC + HOVAL, columbus, listw)
loobj1 <- impacts(lobj1, tr=trMatc, R=200)
summary(loobj1, zstats=TRUE, short=TRUE)
lobj1r <- stsls(CRIME ~ INC + HOVAL, columbus, listw, robust=TRUE)
loobj1r <- impacts(lobj1r, tr=trMatc, R=200)
summary(loobj1r, zstats=TRUE, short=TRUE)
lobjIQ5 <- impacts(lobj, tr=trMatc, R=200, Q=5)
summary(lobjIQ5, zstats=TRUE, short=TRUE)
summary(lobjIQ5, zstats=TRUE, short=TRUE, reportQ=TRUE)
impacts(mobj, listw=listw)
impacts(mobj, tr=trMatc)
impacts(mobj, tr=trMC)
summary(impacts(mobj, tr=trMatc, R=200), zstats=TRUE)
## Not run: 
##D mobj1 <- lagsarlm(CRIME ~ INC + HOVAL, columbus, listw, type="mixed", 
##D method="Matrix", fdHess=TRUE)
##D summary(mobj1)
##D summary(impacts(mobj1, tr=trMatc, R=1000), zstats=TRUE, short=TRUE)
##D summary(impacts(mobj, tr=trMatc, R=1000), zstats=TRUE, short=TRUE)
##D mobj2 <- lagsarlm(CRIME ~ INC + HOVAL, columbus, listw, type="mixed", 
##D method="Matrix", fdHess=TRUE, optimHess=TRUE)
##D summary(impacts(mobj2, tr=trMatc, R=1000), zstats=TRUE, short=TRUE)
##D \dontrun{
##D mobj3 <- lagsarlm(CRIME ~ INC + HOVAL, columbus, listw, type="mixed", 
##D method="spam", fdHess=TRUE)
##D summary(impacts(mobj3, tr=trMatc, R=1000), zstats=TRUE, short=TRUE)
##D }
##D data(boston)
##D Wb <- as(as_dgRMatrix_listw(nb2listw(boston.soi)), "CsparseMatrix")
##D trMatb <- trW(Wb, type="mult")
##D gp2mMi <- lagsarlm(log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + 
##D I(RM^2) +  AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT), 
##D data=boston.c, nb2listw(boston.soi), type="mixed", method="Matrix", 
##D fdHess=TRUE, trs=trMatb)
##D summary(gp2mMi)
##D summary(impacts(gp2mMi, tr=trMatb, R=1000), zstats=TRUE, short=TRUE)
##D data(house)
##D lw <- nb2listw(LO_nb)
##D form <- formula(log(price) ~ age + I(age^2) + I(age^3) + log(lotsize) +
##D    rooms + log(TLA) + beds + syear)
##D lobj <- lagsarlm(form, house, lw, method="Matrix",
##D  fdHess=TRUE, trs=trMat)
##D summary(lobj)
##D loobj <- impacts(lobj, tr=trMat, R=1000)
##D summary(loobj, zstats=TRUE, short=TRUE)
##D lobj1 <- stsls(form, house, lw)
##D loobj1 <- impacts(lobj1, tr=trMat, R=1000)
##D summary(loobj1, zstats=TRUE, short=TRUE)
##D mobj <- lagsarlm(form, house, lw, type="mixed",
##D  method="Matrix", fdHess=TRUE, trs=trMat)
##D summary(mobj)
##D moobj <- impacts(mobj, tr=trMat, R=1000)
##D summary(moobj, zstats=TRUE, short=TRUE)
## End(Not run)



cleanEx()
nameEx("include.self")
### * include.self

flush(stderr()); flush(stdout())

### Name: include.self
### Title: Include self in neighbours list
### Aliases: include.self
### Keywords: spatial

### ** Examples

example(columbus)
coords <- coordinates(columbus)
summary(col.gal.nb, coords)
summary(include.self(col.gal.nb), coords)



cleanEx()
nameEx("invIrM")
### * invIrM

flush(stderr()); flush(stdout())

### Name: invIrM
### Title: Compute SAR generating operator
### Aliases: invIrM invIrW powerWeights
### Keywords: spatial

### ** Examples


nb7rt <- cell2nb(7, 7, torus=TRUE)
set.seed(1)
x <- matrix(rnorm(500*length(nb7rt)), nrow=length(nb7rt))
res0 <- apply(invIrM(nb7rt, rho=0.0, method="chol",
 feasible=TRUE) %*% x, 2, function(x) var(x)/length(x))
res2 <- apply(invIrM(nb7rt, rho=0.2, method="chol",
 feasible=TRUE) %*% x, 2, function(x) var(x)/length(x))
res4 <- apply(invIrM(nb7rt, rho=0.4, method="chol",
 feasible=TRUE) %*% x, 2, function(x) var(x)/length(x))
res6 <- apply(invIrM(nb7rt, rho=0.6, method="chol",
 feasible=TRUE) %*% x, 2, function(x) var(x)/length(x))
res8 <- apply(invIrM(nb7rt, rho=0.8, method="chol",
 feasible=TRUE) %*% x, 2, function(x) var(x)/length(x))
res9 <- apply(invIrM(nb7rt, rho=0.9, method="chol",
 feasible=TRUE) %*% x, 2, function(x) var(x)/length(x))
plot(density(res9), col="red", xlim=c(-0.01, max(density(res9)$x)),
  ylim=range(density(res0)$y),
  xlab="estimated variance of the mean",
  main=expression(paste("Effects of spatial autocorrelation for different ",
    rho, " values")))
lines(density(res0), col="black")
lines(density(res2), col="brown")
lines(density(res4), col="green")
lines(density(res6), col="orange")
lines(density(res8), col="pink")
legend(c(-0.02, 0.01), c(7, 25),
 legend=c("0.0", "0.2", "0.4", "0.6", "0.8", "0.9"),
 col=c("black", "brown", "green", "orange", "pink", "red"), lty=1, bty="n")
## Not run: 
##D x <- matrix(rnorm(length(nb7rt)), ncol=1)
##D system.time(e <- invIrM(nb7rt, rho=0.9, method="chol", feasible=TRUE) %*% x)
##D system.time(e <- invIrM(nb7rt, rho=0.9, method="chol", feasible=NULL) %*% x)
##D system.time(e <- invIrM(nb7rt, rho=0.9, method="solve", feasible=TRUE) %*% x)
##D system.time(e <- invIrM(nb7rt, rho=0.9, method="solve", feasible=NULL) %*% x)
##D W <- as(as_dgRMatrix_listw(nb2listw(nb7rt)), "CsparseMatrix")
##D system.time(ee <- powerWeights(W, rho=0.9, X=x))
##D all.equal(e, as(ee, "matrix"), check.attributes=FALSE)
##D nb60rt <- cell2nb(60, 60, torus=TRUE)
##D W <- as(as_dgRMatrix_listw(nb2listw(nb60rt)), "CsparseMatrix")
##D set.seed(1)
##D x <- matrix(rnorm(dim(W)[1]), ncol=1)
##D system.time(ee <- powerWeights(W, rho=0.3, X=x))
##D str(as(ee, "matrix"))
##D obj <- errorsarlm(as(ee, "matrix")[,1] ~ 1, listw=nb2listw(nb60rt), method="Matrix")
##D coefficients(obj)
## End(Not run)



cleanEx()
nameEx("joincount.mc")
### * joincount.mc

flush(stderr()); flush(stdout())

### Name: joincount.mc
### Title: Permutation test for same colour join count statistics
### Aliases: joincount.mc
### Keywords: spatial

### ** Examples

data(oldcol)
HICRIME <- cut(COL.OLD$CRIME, breaks=c(0,35,80), labels=c("low","high"))
names(HICRIME) <- rownames(COL.OLD)
joincount.mc(HICRIME, nb2listw(COL.nb, style="B"), nsim=99)
joincount.test(HICRIME, nb2listw(COL.nb, style="B"))



cleanEx()
nameEx("joincount.multi")
### * joincount.multi

flush(stderr()); flush(stdout())

### Name: joincount.multi
### Title: BB, BW and Jtot join count statistic for k-coloured factors
### Aliases: joincount.multi print.jcmulti
### Keywords: spatial

### ** Examples

data(oldcol)
HICRIME <- cut(COL.OLD$CRIME, breaks=c(0,35,80), labels=c("low","high"))
names(HICRIME) <- rownames(COL.OLD)
joincount.multi(HICRIME, nb2listw(COL.nb, style="B"))
## Not run: 
##D data(hopkins)
##D image(1:32, 1:32, hopkins[5:36,36:5], breaks=c(-0.5, 3.5, 20),
##D  col=c("white", "black"))
##D box()
##D hopkins.rook.nb <- cell2nb(32, 32, type="rook")
##D unlist(spweights.constants(nb2listw(hopkins.rook.nb, style="B")))
##D hopkins.queen.nb <- cell2nb(32, 32, type="queen")
##D hopkins.bishop.nb <- diffnb(hopkins.rook.nb, hopkins.queen.nb, verbose=FALSE)
##D hopkins4 <- hopkins[5:36,36:5]
##D hopkins4[which(hopkins4 > 3, arr.ind=TRUE)] <- 4
##D hopkins4.f <- factor(hopkins4)
##D table(hopkins4.f)
##D joincount.multi(hopkins4.f, nb2listw(hopkins.rook.nb, style="B"))
##D cat("replicates Upton & Fingleton table 3.4 (p. 166)\n")
##D joincount.multi(hopkins4.f, nb2listw(hopkins.bishop.nb, style="B"))
##D cat("replicates Upton & Fingleton table 3.6 (p. 168)\n")
##D joincount.multi(hopkins4.f, nb2listw(hopkins.queen.nb, style="B"))
##D cat("replicates Upton & Fingleton table 3.7 (p. 169)\n")
## End(Not run)



cleanEx()
nameEx("joincount.test")
### * joincount.test

flush(stderr()); flush(stdout())

### Name: joincount.test
### Title: BB join count statistic for k-coloured factors
### Aliases: joincount.test print.jclist
### Keywords: spatial

### ** Examples

data(oldcol)
HICRIME <- cut(COL.OLD$CRIME, breaks=c(0,35,80), labels=c("low","high"))
names(HICRIME) <- rownames(COL.OLD)
joincount.test(HICRIME, nb2listw(COL.nb, style="B"))
joincount.test(HICRIME, nb2listw(COL.nb, style="C"))
joincount.test(HICRIME, nb2listw(COL.nb, style="S"))
joincount.test(HICRIME, nb2listw(COL.nb, style="W"))
by(card(COL.nb), HICRIME, summary)
print(is.symmetric.nb(COL.nb))
coords.OLD <- cbind(COL.OLD$X, COL.OLD$Y)
COL.k4.nb <- knn2nb(knearneigh(coords.OLD, 4))
print(is.symmetric.nb(COL.k4.nb))
joincount.test(HICRIME, nb2listw(COL.k4.nb, style="B"))
cat("Note non-symmetric weights matrix - use listw2U()\n")
joincount.test(HICRIME, listw2U(nb2listw(COL.k4.nb, style="B")))



cleanEx()
nameEx("knearneigh")
### * knearneigh

flush(stderr()); flush(stdout())

### Name: knearneigh
### Title: K nearest neighbours for spatial weights
### Aliases: knearneigh
### Keywords: spatial

### ** Examples

example(columbus)
coords <- coordinates(columbus)
col.knn <- knearneigh(coords, k=4)
plot(columbus, border="grey")
plot(knn2nb(col.knn), coords, add=TRUE)
title(main="K nearest neighbours, k = 4")
data(state)
us48.fipsno <- read.geoda(system.file("etc/weights/us48.txt",
 package="spdep")[1])
if (as.numeric(paste(version$major, version$minor, sep="")) < 19) {
 m50.48 <- match(us48.fipsno$"State.name", state.name)
} else {
 m50.48 <- match(us48.fipsno$"State_name", state.name)
}
xy <- as.matrix(as.data.frame(state.center))[m50.48,]
llk4.nb <- knn2nb(knearneigh(xy, k=4, longlat=FALSE))
gck4.nb <- knn2nb(knearneigh(xy, k=4, longlat=TRUE))
plot(llk4.nb, xy)
plot(diffnb(llk4.nb, gck4.nb), xy, add=TRUE, col="red", lty=2)
title(main="Differences between Euclidean and Great Circle k=4 neighbours")
summary(llk4.nb, xy, longlat=TRUE)
summary(gck4.nb, xy, longlat=TRUE)

xy1 <- SpatialPoints((as.data.frame(state.center))[m50.48,],
  proj4string=CRS("+proj=longlat"))
gck4a.nb <- knn2nb(knearneigh(xy1, k=4))
summary(gck4a.nb, xy1)



cleanEx()
nameEx("knn2nb")
### * knn2nb

flush(stderr()); flush(stdout())

### Name: knn2nb
### Title: Neighbours list from knn object
### Aliases: knn2nb
### Keywords: spatial

### ** Examples

example(columbus)
coords <- coordinates(columbus)
col.knn <- knearneigh(coords, k=4)
plot(columbus, border="grey")
plot(knn2nb(col.knn), coords, add=TRUE)
title(main="K nearest neighbours, k = 4")



cleanEx()
nameEx("lag.listw")
### * lag.listw

flush(stderr()); flush(stdout())

### Name: lag.listw
### Title: Spatial lag of a numeric vector
### Aliases: lag.listw
### Keywords: spatial

### ** Examples

data(oldcol)
Vx <- lag.listw(nb2listw(COL.nb, style="W"), COL.OLD$CRIME)
plot(Vx, COL.OLD$CRIME)
plot(ecdf(COL.OLD$CRIME))
plot(ecdf(Vx), add=TRUE, col.points="red", col.hor="red")
is.na(COL.OLD$CRIME[5]) <- TRUE
VxNA <- lag.listw(nb2listw(COL.nb, style="W"), COL.OLD$CRIME, NAOK=TRUE)



cleanEx()
nameEx("lagmess")
### * lagmess

flush(stderr()); flush(stdout())

### Name: lagmess
### Title: Matrix exponential spatial lag model
### Aliases: lagmess print.lagmess print.summary.lagmess summary.lagmess
###   residuals.lagmess deviance.lagmess coef.lagmess fitted.lagmess
###   logLik.lagmess
### Keywords: spatial

### ** Examples

data(baltimore)
baltimore$AGE <- ifelse(baltimore$AGE < 1, 1, baltimore$AGE)
lw <- nb2listw(knn2nb(knearneigh(cbind(baltimore$X, baltimore$Y), k=7)))
obj1 <- lm(log(PRICE) ~ PATIO + log(AGE) + log(SQFT) + lag(lw, log(AGE)),
 data=baltimore)
lm.morantest(obj1, lw)
lm.LMtests(obj1, lw, test="all")
obj2 <- lagmess(log(PRICE) ~ PATIO + log(AGE) + log(SQFT) + 
 lag(lw, log(AGE)), data=baltimore, listw=lw)
summary(obj2)
obj3 <- lagsarlm(log(PRICE) ~ PATIO + log(AGE) + log(SQFT) + 
 lag(lw, log(AGE)), data=baltimore, listw=lw)
summary(obj3)
data(boston)
lw <- nb2listw(boston.soi)
gp2 <- lagsarlm(log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + I(RM^2)
 +  AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT),
 data=boston.c, lw, method="Matrix")
summary(gp2)
gp2a <- lagmess(CMEDV ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + I(RM^2)
 +  AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT),
 data=boston.c, lw)
summary(gp2a)



cleanEx()
nameEx("lagsarlm")
### * lagsarlm

flush(stderr()); flush(stdout())

### Name: lagsarlm
### Title: Spatial simultaneous autoregressive lag model estimation
### Aliases: lagsarlm
### Keywords: spatial

### ** Examples

data(oldcol)
COL.lag.eig <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="W"), method="eigen", quiet=FALSE)
summary(COL.lag.eig, correlation=TRUE)
COL.lag.eig$fdHess
COL.lag.eig$resvar
W <- as(as_dgRMatrix_listw(nb2listw(COL.nb)), "CsparseMatrix")
trMatc <- trW(W, type="mult")
COL.lag.eig1 <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="W"), control=list(fdHess=TRUE), trs=trMatc)
COL.lag.eig1$fdHess
system.time(COL.lag.M <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb), method="Matrix", quiet=FALSE))
summary(COL.lag.M)
impacts(COL.lag.M, listw=nb2listw(COL.nb))
## Not run: 
##D system.time(COL.lag.sp <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
##D  nb2listw(COL.nb), method="spam", quiet=FALSE))
##D summary(COL.lag.sp)
## End(Not run)
COL.lag.B <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="B"))
summary(COL.lag.B, correlation=TRUE)
COL.mixed.B <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="B"), type="mixed", tol.solve=1e-9)
summary(COL.mixed.B, correlation=TRUE)
COL.mixed.W <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb, style="W"), type="mixed")
summary(COL.mixed.W, correlation=TRUE)
NA.COL.OLD <- COL.OLD
NA.COL.OLD$CRIME[20:25] <- NA
COL.lag.NA <- lagsarlm(CRIME ~ INC + HOVAL, data=NA.COL.OLD,
 nb2listw(COL.nb), na.action=na.exclude, 
 control=list(tol.opt=.Machine$double.eps^0.4))
COL.lag.NA$na.action
COL.lag.NA
resid(COL.lag.NA)
data(boston)
gp2mM <- lagsarlm(log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + 
I(RM^2) +  AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT), 
data=boston.c, nb2listw(boston.soi), type="mixed", method="Matrix")
summary(gp2mM)
W <- as(as_dgRMatrix_listw(nb2listw(boston.soi)), "CsparseMatrix")
trMatb <- trW(W, type="mult")
gp2mMi <- lagsarlm(log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + 
I(RM^2) +  AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT), 
data=boston.c, nb2listw(boston.soi), type="mixed", method="Matrix", 
trs=trMatb)
summary(gp2mMi)



cleanEx()
nameEx("listw2sn")
### * listw2sn

flush(stderr()); flush(stdout())

### Name: listw2sn
### Title: Spatial neighbour sparse representation
### Aliases: listw2sn sn2listw as.spam.listw
### Keywords: spatial

### ** Examples

example(columbus)
col.listw <- nb2listw(col.gal.nb)
col.listw$neighbours[[1]]
col.listw$weights[[1]]
col.sn <- listw2sn(col.listw)
str(col.sn)
## Not run: 
##D col.sp <- as.spam.listw(col.listw)
##D str(col.sp)
## End(Not run)



cleanEx()
nameEx("lm.LMtests")
### * lm.LMtests

flush(stderr()); flush(stdout())

### Name: lm.LMtests
### Title: Lagrange Multiplier diagnostics for spatial dependence in linear
###   models
### Aliases: lm.LMtests print.LMtestlist
### Keywords: spatial

### ** Examples

data(oldcol)
oldcrime.lm <- lm(CRIME ~ HOVAL + INC, data = COL.OLD)
summary(oldcrime.lm)
lm.LMtests(oldcrime.lm, nb2listw(COL.nb), test=c("LMerr", "LMlag", "RLMerr",
  "RLMlag", "SARMA"))
lm.LMtests(oldcrime.lm, nb2listw(COL.nb))
lm.LMtests(residuals(oldcrime.lm), nb2listw(COL.nb))



cleanEx()
nameEx("lm.morantest")
### * lm.morantest

flush(stderr()); flush(stdout())

### Name: lm.morantest
### Title: Moran's I test for residual spatial autocorrelation
### Aliases: lm.morantest listw2U
### Keywords: spatial

### ** Examples

data(oldcol)
oldcrime1.lm <- lm(CRIME ~ 1, data = COL.OLD)
oldcrime.lm <- lm(CRIME ~ HOVAL + INC, data = COL.OLD)
lm.morantest(oldcrime.lm, nb2listw(COL.nb, style="W"))
lm.LMtests(oldcrime.lm, nb2listw(COL.nb, style="W"))
lm.morantest(oldcrime.lm, nb2listw(COL.nb, style="S"))
lm.morantest(oldcrime1.lm, nb2listw(COL.nb, style="W"))
moran.test(COL.OLD$CRIME, nb2listw(COL.nb, style="W"),
 randomisation=FALSE)
oldcrime.wlm <- lm(CRIME ~ HOVAL + INC, data = COL.OLD,
 weights = I(1/AREA))
lm.morantest(oldcrime.wlm, nb2listw(COL.nb, style="W"),
 resfun=weighted.residuals)
lm.morantest(oldcrime.wlm, nb2listw(COL.nb, style="W"),
 resfun=rstudent)
if (require(boot)) {
  oldcrime.lmx <- lm(CRIME ~ HOVAL + INC, data = COL.OLD, x=TRUE)
  listw <- nb2listw(COL.nb, style="W")
  MoraneI.boot <- function(var, i, ...) {
    var <- var[i]
    lmres <- lm(var ~ oldcrime.lmx$x - 1)
    return(moran(x=residuals(lmres), ...)$I)
  }
  boot1 <- boot(residuals(oldcrime.lmx), statistic=MoraneI.boot, R=499,
    sim="permutation", listw=listw, n=length(listw$neighbours),
    S0=Szero(listw))
  zi <- (boot1$t0 - mean(boot1$t))/sqrt(var(boot1$t))
  boot1
  plot(boot1)
  cat("Bootstrap permutation standard deviate:", zi, "\n\n")
  lm.morantest(oldcrime.lm, nb2listw(COL.nb, style="W"))
}



cleanEx()
nameEx("lm.morantest.exact")
### * lm.morantest.exact

flush(stderr()); flush(stdout())

### Name: lm.morantest.exact
### Title: Exact global Moran's I test
### Aliases: lm.morantest.exact print.moranex
### Keywords: spatial

### ** Examples

eire <- readShapePoly(system.file("etc/shapes/eire.shp", package="spdep")[1],
  ID="names", proj4string=CRS("+proj=utm +zone=30 +units=km"))
eire.nb <- poly2nb(eire)
#data(eire)
e.lm <- lm(OWNCONS ~ ROADACC, data=eire)
lm.morantest(e.lm, nb2listw(eire.nb))
lm.morantest.sad(e.lm, nb2listw(eire.nb))
lm.morantest.exact(e.lm, nb2listw(eire.nb))
lm.morantest.exact(e.lm, nb2listw(eire.nb), useTP=TRUE)



cleanEx()
nameEx("lm.morantest.sad")
### * lm.morantest.sad

flush(stderr()); flush(stdout())

### Name: lm.morantest.sad
### Title: Saddlepoint approximation of global Moran's I test
### Aliases: lm.morantest.sad print.moransad summary.moransad
###   print.summary.moransad
### Keywords: spatial

### ** Examples

eire <- readShapePoly(system.file("etc/shapes/eire.shp", package="spdep")[1],
  ID="names", proj4string=CRS("+proj=utm +zone=30 +units=km"))
eire.nb <- poly2nb(eire)
#data(eire)
e.lm <- lm(OWNCONS ~ ROADACC, data=eire)
lm.morantest(e.lm, nb2listw(eire.nb))
lm.morantest.sad(e.lm, nb2listw(eire.nb))
summary(lm.morantest.sad(e.lm, nb2listw(eire.nb)))
e.wlm <- lm(OWNCONS ~ ROADACC, data=eire, weights=RETSALE)
lm.morantest(e.wlm, nb2listw(eire.nb), resfun=rstudent)
lm.morantest.sad(e.wlm, nb2listw(eire.nb), resfun=rstudent)



cleanEx()
nameEx("localG")
### * localG

flush(stderr()); flush(stdout())

### Name: localG
### Title: G and Gstar local spatial statistics
### Aliases: localG
### Keywords: spatial

### ** Examples

data(getisord)
xycoords <- cbind(xyz$x, xyz$y)
nb30 <- dnearneigh(xycoords, 0, 30)
G30 <- localG(xyz$val, nb2listw(nb30, style="B"))
G30[length(xyz$val)-136]
nb60 <- dnearneigh(xycoords, 0, 60)
G60 <- localG(xyz$val, nb2listw(nb60, style="B"))
G60[length(xyz$val)-136]
nb90 <- dnearneigh(xycoords, 0, 90)
G90 <- localG(xyz$val, nb2listw(nb90, style="B"))
G90[length(xyz$val)-136]
nb120 <- dnearneigh(xycoords, 0, 120)
G120 <- localG(xyz$val, nb2listw(nb120, style="B"))
G120[length(xyz$val)-136]
nb150 <- dnearneigh(xycoords, 0, 150)
G150 <- localG(xyz$val, nb2listw(nb150, style="B"))
G150[length(xyz$val)-136]
brks <- seq(-5,5,1)
cm.col <- cm.colors(length(brks)-1)
image(x, y, t(matrix(G30, nrow=16, ncol=16, byrow=TRUE)),
  breaks=brks, col=cm.col, asp=1)
text(xyz$x, xyz$y, round(G30, digits=1), cex=0.7)
polygon(c(195,225,225,195), c(195,195,225,225), lwd=2)
title(main=expression(paste("Values of the ", G[i], " statistic")))
G30s <- localG(xyz$val, nb2listw(include.self(nb30),
 style="B"))
cat("value according to Getis and Ord's eq. 14.2, p. 263 (1996)\n")
G30s[length(xyz$val)-136]
cat(paste("value given by Getis and Ord (1996), p. 267",
  "(division by n-1 rather than n \n in variance)\n"))
G30s[length(xyz$val)-136] *
  (sqrt(sum(scale(xyz$val, scale=FALSE)^2)/length(xyz$val)) /
  sqrt(var(xyz$val)))
image(x, y, t(matrix(G30s, nrow=16, ncol=16, byrow=TRUE)),
  breaks=brks, col=cm.col, asp=1)
text(xyz$x, xyz$y, round(G30s, digits=1), cex=0.7)
polygon(c(195,225,225,195), c(195,195,225,225), lwd=2)
title(main=expression(paste("Values of the ", G[i]^"*", " statistic")))



cleanEx()
nameEx("localmoran")
### * localmoran

flush(stderr()); flush(stdout())

### Name: localmoran
### Title: Local Moran's I statistic
### Aliases: localmoran
### Keywords: spatial

### ** Examples

data(afcon)
oid <- order(afcon$id)
resI <- localmoran(afcon$totcon, nb2listw(paper.nb))
printCoefmat(data.frame(resI[oid,], row.names=afcon$name[oid]),
 check.names=FALSE)
hist(resI[,5])
resI <- localmoran(afcon$totcon, nb2listw(paper.nb),
 p.adjust.method="bonferroni")
printCoefmat(data.frame(resI[oid,], row.names=afcon$name[oid]),
 check.names=FALSE)
hist(resI[,5])
totcon <-afcon$totcon
is.na(totcon) <- sample(1:length(totcon), 5)
totcon
resI.na <- localmoran(totcon, nb2listw(paper.nb), na.action=na.exclude,
 zero.policy=TRUE)
if (class(attr(resI.na, "na.action")) == "exclude") {
 print(data.frame(resI.na[oid,], row.names=afcon$name[oid]), digits=2)
} else print(resI.na, digits=2)
resG <- localG(afcon$totcon, nb2listw(include.self(paper.nb)))
print(data.frame(resG[oid], row.names=afcon$name[oid]), digits=2)




cleanEx()
nameEx("localmoran.exact")
### * localmoran.exact

flush(stderr()); flush(stdout())

### Name: localmoran.exact
### Title: Exact local Moran's Ii tests
### Aliases: localmoran.exact localmoran.exact.alt print.localmoranex
###   as.data.frame.localmoranex
### Keywords: spatial

### ** Examples

eire <- readShapePoly(system.file("etc/shapes/eire.shp", package="spdep")[1],
  ID="names", proj4string=CRS("+proj=utm +zone=30 +units=km"))
eire.nb <- poly2nb(eire)
#data(eire)
e.lm <- lm(OWNCONS ~ ROADACC, data=eire)
localmoran.sad(e.lm, nb=eire.nb)
localmoran.exact(e.lm, nb=eire.nb)
localmoran.exact(e.lm, nb=eire.nb, useTP=TRUE)
e.errorsar <- errorsarlm(OWNCONS ~ ROADACC, data=eire,
 listw=nb2listw(eire.nb))
lm.target <- lm(e.errorsar$tary ~ e.errorsar$tarX - 1)
localmoran.exact.alt(lm.target, nb=eire.nb)
Omega <- invIrW(nb2listw(eire.nb), rho=0.6)
Omega1 <- tcrossprod(Omega)
localmoran.exact.alt(lm.target, nb=eire.nb, Omega=Omega1)
localmoran.exact.alt(lm.target, nb=eire.nb, Omega=Omega1, useTP=TRUE)



cleanEx()
nameEx("localmoran.sad")
### * localmoran.sad

flush(stderr()); flush(stdout())

### Name: localmoran.sad
### Title: Saddlepoint approximation of local Moran's Ii tests
### Aliases: localmoran.sad listw2star print.summary.localmoransad
###   summary.localmoransad print.localmoransad as.data.frame.localmoransad
### Keywords: spatial

### ** Examples

eire <- readShapePoly(system.file("etc/shapes/eire.shp", package="spdep")[1],
  ID="names", proj4string=CRS("+proj=utm +zone=30 +units=km"))
eire.nb <- poly2nb(eire)
#data(eire)
e.lm <- lm(OWNCONS ~ ROADACC, data=eire)
e.locmor <- summary(localmoran.sad(e.lm, nb=eire.nb))
e.locmor
mean(e.locmor[,1])
lm.morantest(e.lm, nb2listw(eire.nb))
hist(e.locmor[,"Pr. (Sad)"])
e.wlm <- lm(OWNCONS ~ ROADACC, data=eire, weights=RETSALE)
e.locmorw1 <- summary(localmoran.sad(e.wlm, nb=eire.nb, resfun=weighted.residuals))
e.locmorw1
e.locmorw2 <- summary(localmoran.sad(e.wlm, nb=eire.nb, resfun=rstudent))
e.locmorw2
e.errorsar <- errorsarlm(OWNCONS ~ ROADACC, data=eire,
  listw=nb2listw(eire.nb))
e.errorsar
lm.target <- lm(e.errorsar$tary ~ e.errorsar$tarX - 1)
e.clocmor <- summary(localmoran.sad(lm.target, nb=eire.nb))
e.clocmor
hist(e.clocmor[,"Pr. (Sad)"])



cleanEx()
nameEx("mat2listw")
### * mat2listw

flush(stderr()); flush(stdout())

### Name: mat2listw
### Title: Convert a square spatial weights matrix to a weights list object
### Aliases: mat2listw
### Keywords: spatial

### ** Examples

example(columbus)
coords <- coordinates(columbus)
col005 <- dnearneigh(coords, 0, 0.5, attr(col.gal.nb, "region.id"))
summary(col005)
col005.w.mat <- nb2mat(col005, zero.policy=TRUE)
col005.w.b <- mat2listw(col005.w.mat)
summary(col005.w.b$neighbours)
diffnb(col005, col005.w.b$neighbours)
col005.w.mat.3T <- kronecker(diag(3), col005.w.mat)
col005.w.b.3T <- mat2listw(col005.w.mat.3T, style="W")
summary(col005.w.b.3T$neighbours)
W <- as(as_dgRMatrix_listw(nb2listw(col005, style="W", zero.policy=TRUE)), "CsparseMatrix")
col005.spM <- mat2listw(W)
summary(col005.spM$neighbours)
diffnb(col005, col005.spM$neighbours)
IW <- kronecker(Diagonal(3), W)
col005.spM.3T <- mat2listw(IW, style="W")
summary(col005.spM.3T$neighbours)



cleanEx()
nameEx("moran")
### * moran

flush(stderr()); flush(stdout())

### Name: moran
### Title: Compute Moran's I
### Aliases: moran
### Keywords: spatial

### ** Examples

data(oldcol)
col.W <- nb2listw(COL.nb, style="W")
crime <- COL.OLD$CRIME
str(moran(crime, col.W, length(COL.nb), Szero(col.W)))
is.na(crime) <- sample(1:length(crime), 10)
str(moran(crime, col.W, length(COL.nb), Szero(col.W), NAOK=TRUE))



cleanEx()
nameEx("moran.mc")
### * moran.mc

flush(stderr()); flush(stdout())

### Name: moran.mc
### Title: Permutation test for Moran's I statistic
### Aliases: moran.mc
### Keywords: spatial

### ** Examples

data(oldcol)
colw <- nb2listw(COL.nb, style="W")
nsim <- 99
set.seed(1234)
sim1 <- moran.mc(COL.OLD$CRIME, listw=colw, nsim=nsim)
sim1
mean(sim1$res[1:nsim])
var(sim1$res[1:nsim])
summary(sim1$res[1:nsim])
MoranI.boot <- function(var, i, ...) {
var <- var[i]
return(moran(x=var, ...)$I)
}
set.seed(1234)
library(boot)
boot1 <- boot(COL.OLD$CRIME, statistic=MoranI.boot, R=nsim,
 sim="permutation", listw=colw, n=nrow(COL.OLD), S0=Szero(colw))
boot1
plot(boot1)
mean(boot1$t)
var(boot1$t)
summary(boot1$t)
colold.lags <- nblag(COL.nb, 3)
set.seed(1234)
sim2 <- moran.mc(COL.OLD$CRIME, nb2listw(colold.lags[[2]],
 style="W"), nsim=nsim)
summary(sim2$res[1:nsim])
sim3 <- moran.mc(COL.OLD$CRIME, nb2listw(colold.lags[[3]],
 style="W"), nsim=nsim)
summary(sim3$res[1:nsim])



cleanEx()
nameEx("moran.plot")
### * moran.plot

flush(stderr()); flush(stdout())

### Name: moran.plot
### Title: Moran scatterplot
### Aliases: moran.plot
### Keywords: spatial

### ** Examples

data(afcon)
moran.plot(afcon$totcon, nb2listw(paper.nb),
 labels=as.character(afcon$name), pch=19)
moran.plot(as.vector(scale(afcon$totcon)), nb2listw(paper.nb),
 labels=as.character(afcon$name), xlim=c(-2, 4), ylim=c(-2,4), pch=19)



cleanEx()
nameEx("moran.test")
### * moran.test

flush(stderr()); flush(stdout())

### Name: moran.test
### Title: Moran's I test for spatial autocorrelation
### Aliases: moran.test
### Keywords: spatial

### ** Examples

data(oldcol)
coords.OLD <- cbind(COL.OLD$X, COL.OLD$Y)
moran.test(COL.OLD$CRIME, nb2listw(COL.nb, style="W"))
moran.test(COL.OLD$CRIME, nb2listw(COL.nb, style="B"))
moran.test(COL.OLD$CRIME, nb2listw(COL.nb, style="C"))
moran.test(COL.OLD$CRIME, nb2listw(COL.nb, style="S"))
moran.test(COL.OLD$CRIME, nb2listw(COL.nb, style="W"),
 randomisation=FALSE)
colold.lags <- nblag(COL.nb, 3)
moran.test(COL.OLD$CRIME, nb2listw(colold.lags[[2]],
 style="W"))
moran.test(COL.OLD$CRIME, nb2listw(colold.lags[[3]],
 style="W"))
print(is.symmetric.nb(COL.nb))
COL.k4.nb <- knn2nb(knearneigh(coords.OLD, 4))
print(is.symmetric.nb(COL.k4.nb))
moran.test(COL.OLD$CRIME, nb2listw(COL.k4.nb, style="W"))
moran.test(COL.OLD$CRIME, nb2listw(COL.k4.nb, style="W"),
 randomisation=FALSE)
cat("Note: non-symmetric weights matrix, use listw2U()")
moran.test(COL.OLD$CRIME, listw2U(nb2listw(COL.k4.nb,
 style="W")))
moran.test(COL.OLD$CRIME, listw2U(nb2listw(COL.k4.nb,
 style="W")), randomisation=FALSE)
ranks <- rank(COL.OLD$CRIME)
names(ranks) <- rownames(COL.OLD)
moran.test(ranks, nb2listw(COL.nb, style="W"), rank=TRUE)
crime <- COL.OLD$CRIME
is.na(crime) <- sample(1:length(crime), 10)
res <- try(moran.test(crime, nb2listw(COL.nb, style="W"),
 na.action=na.fail))
res
moran.test(crime, nb2listw(COL.nb, style="W"), zero.policy=TRUE,
 na.action=na.omit)
moran.test(crime, nb2listw(COL.nb, style="W"), zero.policy=TRUE,
 na.action=na.exclude)
moran.test(crime, nb2listw(COL.nb, style="W"), na.action=na.pass)



cleanEx()
nameEx("mstree")
### * mstree

flush(stderr()); flush(stdout())

### Name: mstree
### Title: Find the minimal spanning tree
### Aliases: mstree
### Keywords: graphs spatial

### ** Examples

### loading data
bh <- readShapePoly(system.file("etc/shapes/bhicv.shp",
      package="spdep")[1])
### data padronized
dpad <- data.frame(scale(bh@data[,5:8]))

### neighboorhod list 
bh.nb <- poly2nb(bh)

### calculing costs
lcosts <- nbcosts(bh.nb, dpad)

### making listw
nb.w <- nb2listw(bh.nb, lcosts, style="B")

### find a minimum spanning tree
system.time(mst.bh <- mstree(nb.w,5))

dim(mst.bh)

head(mst.bh)
tail(mst.bh)

### the mstree plot
par(mar=c(0,0,0,0))
plot(mst.bh, coordinates(bh), col=2, 
     cex.lab=.7, cex.circles=0.035, fg="blue")
plot(bh, border=gray(.5), add=TRUE)




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("nb2INLA")
### * nb2INLA

flush(stderr()); flush(stdout())

### Name: nb2INLA
### Title: Output spatial neighbours for INLA
### Aliases: nb2INLA
### Keywords: spatial

### ** Examples

example(columbus)
td <- tempdir()
x <- nb2INLA(paste(td, "columbus-INLA.adj", sep="/"), col.gal.nb)



cleanEx()
nameEx("nb2WB")
### * nb2WB

flush(stderr()); flush(stdout())

### Name: nb2WB
### Title: Output spatial weights for WinBUGS
### Aliases: nb2WB listw2WB
### Keywords: spatial

### ** Examples

example(columbus)
x <- nb2WB(col.gal.nb)
dput(x, control=NULL)
x <- listw2WB(nb2listw(col.gal.nb))
dput(x, control=NULL)



cleanEx()
nameEx("nb2blocknb")
### * nb2blocknb

flush(stderr()); flush(stdout())

### Name: nb2blocknb
### Title: Block up neighbour list for location-less observations
### Aliases: nb2blocknb
### Keywords: spatial

### ** Examples

## Not run: 
##D data(boston)
##D summary(as.vector(table(boston.c$TOWN)))
##D townaggr <- aggregate(boston.utm, list(town=boston.c$TOWN), mean)
##D block.rel <- graph2nb(relativeneigh(as.matrix(townaggr[,2:3])),
##D  as.character(townaggr[,1]), sym=TRUE)
##D block.rel
##D print(is.symmetric.nb(block.rel))
##D plot(block.rel, as.matrix(townaggr[,2:3]))
##D points(boston.utm, pch=18, col="lightgreen")
##D block.nb <- nb2blocknb(block.rel, as.character(boston.c$TOWN))
##D block.nb
##D print(is.symmetric.nb(block.nb))
##D plot(block.nb, boston.utm)
##D points(boston.utm, pch=18, col="lightgreen")
##D moran.test(boston.c$CMEDV, nb2listw(boston.soi))
##D moran.test(boston.c$CMEDV, nb2listw(block.nb))
## End(Not run)



cleanEx()
nameEx("nb2lines")
### * nb2lines

flush(stderr()); flush(stdout())

### Name: nb2lines
### Title: Use arc-type shapefiles for import and export of weights
### Aliases: nb2lines listw2lines df2sn
### Keywords: spatial

### ** Examples

example(columbus)
coords <- coordinates(columbus)
res <- listw2lines(nb2listw(col.gal.nb), coords)
summary(res)
fn <- paste(tempdir(), "nbshape", sep="/")
writeLinesShape(res, fn=fn)
inMap <- readShapeLines(fn)
summary(inMap)
diffnb(sn2listw(df2sn(as(inMap, "data.frame")))$neighbours, col.gal.nb)



cleanEx()
nameEx("nb2listw")
### * nb2listw

flush(stderr()); flush(stdout())

### Name: nb2listw
### Title: Spatial weights for neighbours lists
### Aliases: nb2listw
### Keywords: spatial

### ** Examples

example(columbus)
coords <- coordinates(columbus)
cards <- card(col.gal.nb)
col.w <- nb2listw(col.gal.nb)
plot(cards, unlist(lapply(col.w$weights, sum)),xlim=c(0,10),
ylim=c(0,10), xlab="number of links", ylab="row sums of weights")
col.b <- nb2listw(col.gal.nb, style="B")
points(cards, unlist(lapply(col.b$weights, sum)), col="red")
col.c <- nb2listw(col.gal.nb, style="C")
points(cards, unlist(lapply(col.c$weights, sum)), col="green")
col.u <- nb2listw(col.gal.nb, style="U")
points(cards, unlist(lapply(col.u$weights, sum)), col="orange")
col.s <- nb2listw(col.gal.nb, style="S")
points(cards, unlist(lapply(col.s$weights, sum)), col="blue")
legend(x=c(0, 1), y=c(7, 9), legend=c("W", "B", "C", "U", "S"),
col=c("black", "red", "green", "orange", "blue"), pch=rep(1,5))
summary(nb2listw(col.gal.nb, style="minmax"))
dlist <- nbdists(col.gal.nb, coords)
dlist <- lapply(dlist, function(x) 1/x)
col.w.d <- nb2listw(col.gal.nb, glist=dlist)
summary(unlist(col.w$weights))
summary(unlist(col.w.d$weights))
# introducing other conditions into weights - only earlier sales count
# see http://sal.uiuc.edu/pipermail/openspace/2005-October/000610.html
data(baltimore)
set.seed(211)
dates <- sample(1:500, nrow(baltimore), replace=TRUE)
nb_15nn <- knn2nb(knearneigh(cbind(baltimore$X, baltimore$Y), k=15))
glist <- vector(mode="list", length=length(nb_15nn))
for (i in seq(along=nb_15nn))
  glist[[i]] <- ifelse(dates[i] > dates[nb_15nn[[i]]], 1, 0)
listw_15nn_dates <- nb2listw(nb_15nn, glist=glist, style="B")
which(lag(listw_15nn_dates, baltimore$PRICE) == 0.0)
which(sapply(glist, sum) == 0)
ex <- which(sapply(glist, sum) == 0)[1]
dates[ex]
dates[nb_15nn[[ex]]]



cleanEx()
nameEx("nb2mat")
### * nb2mat

flush(stderr()); flush(stdout())

### Name: nb2mat
### Title: Spatial weights matrices for neighbours lists
### Aliases: nb2mat listw2mat
### Keywords: spatial

### ** Examples

example(columbus)
coords <- coordinates(columbus)
col005 <- dnearneigh(coords, 0, 0.5, attr(col.gal.nb, "region.id"))
summary(col005)
col005.w.mat <- nb2mat(col005, zero.policy=TRUE)
table(round(apply(col005.w.mat, 1, sum)))



cleanEx()
nameEx("nbdists")
### * nbdists

flush(stderr()); flush(stdout())

### Name: nbdists
### Title: Spatial link distance measures
### Aliases: nbdists
### Keywords: spatial

### ** Examples

example(columbus)
coords <- coordinates(columbus)
dlist <- nbdists(col.gal.nb, coords)
dlist <- lapply(dlist, function(x) 1/x)
stem(unlist(dlist))



cleanEx()
nameEx("nblag")
### * nblag

flush(stderr()); flush(stdout())

### Name: nblag
### Title: Higher order neighbours lists
### Aliases: nblag nblag_cumul
### Keywords: spatial

### ** Examples

example(columbus)
coords <- coordinates(columbus)
summary(col.gal.nb, coords)
plot(columbus, border="grey")
plot(col.gal.nb, coords, add=TRUE)
title(main="GAL order 1 (black) and 2 (red) links")
col.lags <- nblag(col.gal.nb, 2)
lapply(col.lags, print)
summary(col.lags[[2]], coords)
plot(col.lags[[2]], coords, add=TRUE, col="red", lty=2)
cuml <- nblag_cumul(col.lags)
cuml



cleanEx()
nameEx("nboperations")
### * nboperations

flush(stderr()); flush(stdout())

### Name: nb.set.operations
### Title: Set operations on neighborhood objects
### Aliases: intersect.nb union.nb setdiff.nb complement.nb
### Keywords: spatial

### ** Examples

example(columbus)
coords <- coordinates(columbus)
col.tri.nb <- tri2nb(coords)
oldpar <- par(mfrow=c(1,2))
col.soi.nb <- graph2nb(soi.graph(col.tri.nb, coords))
plot(columbus, border="grey")
plot(col.soi.nb, coords, add=TRUE)
title(main="Sphere of Influence Graph")
plot(columbus, border="grey")
plot(complement.nb(col.soi.nb), coords, add=TRUE)
title(main="Complement of Sphere of Influence Graph")
par(mfrow=c(2,2))
col2 <- droplinks(col.gal.nb, 21)
plot(intersect.nb(col.gal.nb, col2), coords)
title(main="Intersect")
plot(union.nb(col.gal.nb, col2), coords)
title(main="Union")
plot(setdiff.nb(col.gal.nb, col2), coords)
title(main="Set diff")
par(oldpar)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("nc.sids")
### * nc.sids

flush(stderr()); flush(stdout())

### Name: nc.sids
### Title: North Carolina SIDS data
### Aliases: nc.sids ncCR85.nb ncCC89.nb sidspolys sidscents
### Keywords: datasets

### ** Examples

nc.sids <- readShapePoly(system.file("etc/shapes/sids.shp", package="spdep")[1],
  ID="FIPSNO", proj4string=CRS("+proj=longlat +ellps=clrk66"))
rn <- sapply(slot(nc.sids, "polygons"), function(x) slot(x, "ID"))
ncCC89_nb <- read.gal(system.file("etc/weights/ncCC89.gal", package="spdep")[1],
  region.id=rn)
ncCR85_nb <- read.gal(system.file("etc/weights/ncCR85.gal", package="spdep")[1],
  region.id=rn)
## Not run: 
##D plot(nc.sids, border="grey")
##D plot(ncCR85_nb, coordinates(nc.sids), add=TRUE, col="blue")
##D plot(nc.sids, border="grey")
##D plot(ncCC89_nb, coordinates(nc.sids), add=TRUE, col="blue")
## End(Not run)



cleanEx()
nameEx("p.adjustSP")
### * p.adjustSP

flush(stderr()); flush(stdout())

### Name: p.adjustSP
### Title: Adjust local association measures' p-values
### Aliases: p.adjustSP
### Keywords: spatial

### ** Examples

data(afcon)
oid <- order(afcon$id)
resG <- as.vector(localG(afcon$totcon, nb2listw(include.self(paper.nb))))
non <- format.pval(pnorm(2*(abs(resG)), lower.tail=FALSE), 2)
bon <- format.pval(p.adjustSP(pnorm(2*(abs(resG)), lower.tail=FALSE),
 paper.nb, "bonferroni"), 2)
tot <- format.pval(p.adjust(pnorm(2*(abs(resG)), lower.tail=FALSE),
 "bonferroni", n=length(resG)), 2)
data.frame(resG, non, bon, tot, row.names=afcon$name)[oid,]



cleanEx()
nameEx("plot.mst")
### * plot.mst

flush(stderr()); flush(stdout())

### Name: plot.mst
### Title: Plot the Minimum Spanning Tree
### Aliases: plot.mst
### Keywords: hplot tree

### ** Examples

### see example in mstree function documentation



cleanEx()
nameEx("plot.nb")
### * plot.nb

flush(stderr()); flush(stdout())

### Name: plot.nb
### Title: Plot a neighbours list
### Aliases: plot.nb plot.listw
### Keywords: spatial

### ** Examples

example(columbus)
coords <- coordinates(columbus)
plot(col.gal.nb, coords)
title(main="GAL order 1 links with first nearest neighbours in red")
col.knn <- knearneigh(coords, k=1)
plot(knn2nb(col.knn), coords, add=TRUE, col="red", length=0.08)



cleanEx()
nameEx("plot.skater")
### * plot.skater

flush(stderr()); flush(stdout())

### Name: plot.skater
### Title: Plot the object of skater class
### Aliases: plot.skater
### Keywords: hplot cluster

### ** Examples

### see example in the skater function documentation



cleanEx()
nameEx("poly2nb")
### * poly2nb

flush(stderr()); flush(stdout())

### Name: poly2nb
### Title: Construct neighbours list from polygon list
### Aliases: poly2nb
### Keywords: spatial

### ** Examples

example(columbus)
coords <- coordinates(columbus)
xx <- poly2nb(columbus)
dxx <- diffnb(xx, col.gal.nb)
plot(columbus, border="grey")
plot(col.gal.nb, coords, add=TRUE)
plot(dxx, coords, add=TRUE, col="red")
title(main=paste("Differences (red) in Columbus GAL weights (black)",
 "and polygon generated queen weights", sep="\n"))
xxx <- poly2nb(columbus, queen=FALSE)
dxxx <- diffnb(xxx, col.gal.nb)
plot(columbus, border = "grey")
plot(col.gal.nb, coords, add = TRUE)
plot(dxxx, coords, add = TRUE, col = "red")
title(main=paste("Differences (red) in Columbus GAL weights (black)",
 "and polygon generated rook weights", sep="\n"))
cards <- card(xx)
maxconts <- which(cards == max(cards))
if(length(maxconts) > 1) maxconts <- maxconts[1]
fg <- rep("grey", length(cards))
fg[maxconts] <- "red"
fg[xx[[maxconts]]] <- "green"
plot(columbus, col=fg)
title(main="Region with largest number of contiguities")
example(nc.sids)
system.time(xxnb <- poly2nb(nc.sids))
plot(nc.sids)
plot(xxnb, coordinates(nc.sids), add=TRUE, col="blue")



cleanEx()
nameEx("predict.sarlm")
### * predict.sarlm

flush(stderr()); flush(stdout())

### Name: predict.sarlm
### Title: Prediction for spatial simultaneous autoregressive linear model
###   objects
### Aliases: predict.sarlm print.sarlm.pred
### Keywords: spatial

### ** Examples

data(oldcol)
COL.lag.eig <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, nb2listw(COL.nb))
COL.mix.eig <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, nb2listw(COL.nb),
  type="mixed")
COL.err.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, nb2listw(COL.nb))
print(p1 <- predict(COL.mix.eig))
print(p2 <- predict(COL.mix.eig, newdata=COL.OLD, listw=nb2listw(COL.nb)))
AIC(COL.mix.eig)
sqrt(deviance(COL.mix.eig)/length(COL.nb))
sqrt(sum((COL.OLD$CRIME - as.vector(p1))^2)/length(COL.nb))
sqrt(sum((COL.OLD$CRIME - as.vector(p2))^2)/length(COL.nb))
AIC(COL.err.eig)
sqrt(deviance(COL.err.eig)/length(COL.nb))
sqrt(sum((COL.OLD$CRIME - as.vector(predict(COL.err.eig)))^2)/length(COL.nb))
sqrt(sum((COL.OLD$CRIME - as.vector(predict(COL.err.eig, newdata=COL.OLD,
  listw=nb2listw(COL.nb))))^2)/length(COL.nb))
AIC(COL.lag.eig)
sqrt(deviance(COL.lag.eig)/length(COL.nb))
sqrt(sum((COL.OLD$CRIME - as.vector(predict(COL.lag.eig)))^2)/length(COL.nb))
sqrt(sum((COL.OLD$CRIME - as.vector(predict(COL.lag.eig, newdata=COL.OLD,
  listw=nb2listw(COL.nb))))^2)/length(COL.nb))



cleanEx()
nameEx("probmap")
### * probmap

flush(stderr()); flush(stdout())

### Name: probmap
### Title: Probability mapping for rates
### Aliases: probmap
### Keywords: spatial

### ** Examples

example(auckland)
res <- probmap(auckland$M77_85, 9*auckland$Und5_81)
rt <- sum(auckland$M77_85)/sum(9*auckland$Und5_81)
ppois_pmap <- numeric(length(auckland$Und5_81))
for (i in seq(along=ppois_pmap)) {
ppois_pmap[i] <- poisson.test(auckland$M77_85[i], r=rt,
  T=(9*auckland$Und5_81[i]), alternative="less")$p.value
}
all.equal(ppois_pmap, res$pmap)
brks <- c(-Inf,2,2.5,3,3.5,Inf)
cols <- grey(6:2/7)
plot(auckland, col=cols[findInterval(res$raw*1000, brks, all.inside=TRUE)])
legend("bottomleft", fill=cols, legend=leglabs(brks), bty="n")
title(main="Crude (raw) estimates of infant mortality per 1000 per year")
brks <- c(-Inf,47,83,118,154,190,Inf)
cols <- cm.colors(6)
plot(auckland, col=cols[findInterval(res$relRisk, brks, all.inside=TRUE)])
legend("bottomleft", fill=cols, legend=leglabs(brks), bty="n")
title(main="Standardised mortality ratios for Auckland child deaths")
brks <- c(0,0.05,0.1,0.2,0.8,0.9,0.95,1)
cols <- cm.colors(7)
plot(auckland, col=cols[findInterval(res$pmap, brks, all.inside=TRUE)])
legend("bottomleft", fill=cols, legend=leglabs(brks), bty="n")
title(main="Poisson probabilities for Auckland child mortality")



cleanEx()
nameEx("prunecost")
### * prunecost

flush(stderr()); flush(stdout())

### Name: prunecost
### Title: Compute cost of prune each edge
### Aliases: prunecost
### Keywords: graphs cluster

### ** Examples

d <- data.frame(a=-2:2, b=runif(5))
e <- matrix(c(1,2, 2,3, 3,4, 4,5), ncol=2, byrow=TRUE)

sum(sweep(d, 2, colMeans(d))^2)

prunecost(e, d)



cleanEx()
nameEx("prunemst")
### * prunemst

flush(stderr()); flush(stdout())

### Name: prunemst
### Title: Prune a Minimun Spanning Tree
### Aliases: prunemst
### Keywords: tree cluster

### ** Examples

e <- matrix(c(2,3, 1,2, 3,4, 4,5), ncol=2, byrow=TRUE)
e
prunemst(e)
prunemst(e, only.nodes=FALSE)



cleanEx()
nameEx("read.gal")
### * read.gal

flush(stderr()); flush(stdout())

### Name: read.gal
### Title: Read a GAL lattice file into a neighbours list
### Aliases: read.gal read.geoda
### Keywords: spatial

### ** Examples

us48.fipsno <- read.geoda(system.file("etc/weights/us48.txt",
 package="spdep")[1])
us48.q <- read.gal(system.file("etc/weights/us48_q.GAL", package="spdep")[1],
 us48.fipsno$Fipsno)
us48.r <- read.gal(system.file("etc/weights/us48_rk.GAL", package="spdep")[1],
 us48.fipsno$Fipsno)
data(state)
if (as.numeric(paste(version$major, version$minor, sep="")) < 19) {
 m50.48 <- match(us48.fipsno$"State.name", state.name)
} else {
 m50.48 <- match(us48.fipsno$"State_name", state.name)
}
plot(us48.q, as.matrix(as.data.frame(state.center))[m50.48,])
plot(diffnb(us48.r, us48.q),
 as.matrix(as.data.frame(state.center))[m50.48,], add=TRUE, col="red")
title(main="Differences between rook and queen criteria imported neighbours lists")



cleanEx()
nameEx("read.gwt2nb")
### * read.gwt2nb

flush(stderr()); flush(stdout())

### Name: read.gwt2nb
### Title: Read and write spatial neighbour files
### Aliases: read.gwt2nb write.sn2gwt read.dat2listw write.sn2dat
### Keywords: spatial

### ** Examples

data(baltimore)
STATION <- baltimore$STATION
gwt1 <- read.gwt2nb(system.file("etc/weights/baltk4.GWT", package="spdep")[1],
 STATION)
cat(paste("Neighbours list symmetry;", is.symmetric.nb(gwt1, FALSE, TRUE),
 "\n"))
listw1 <- nb2listw(gwt1, style="B", glist=attr(gwt1, "GeoDa")$dist)
tmpGWT <- tempfile()
write.sn2gwt(listw2sn(listw1), tmpGWT)
gwt2 <- read.gwt2nb(tmpGWT, STATION)
cat(paste("Neighbours list symmetry;", is.symmetric.nb(gwt2, FALSE, TRUE),
 "\n"))
diffnb(gwt1, gwt2)
data(oldcol)
tmpMAT <- tempfile()
COL.W <- nb2listw(COL.nb)
write.sn2dat(listw2sn(COL.W), tmpMAT)
listwmat1 <- read.dat2listw(tmpMAT)
diffnb(listwmat1$neighbours, COL.nb, verbose=TRUE)
listwmat2 <- read.dat2listw(system.file("etc/weights/wmat.dat", 
 package="spdep")[1])
diffnb(listwmat1$neighbours, listwmat2$neighbours, verbose=TRUE)



cleanEx()
nameEx("rotation")
### * rotation

flush(stderr()); flush(stdout())

### Name: Rotation
### Title: Rotate a set of point by a certain angle
### Aliases: Rotation
### Keywords: manip

### ** Examples

set.seed(1)
### Create a set of coordinates
coords<-cbind(runif(20),runif(20))

### Create a series of angles
rad<-seq(0,pi,l=20)

opar <- par(mfrow=c(5,4))
for(i in rad){
	coords.rot<-Rotation(coords,i)
	plot(coords.rot)
}
par(opar)

### Rotate the coordinates by an angle of 90 degrees
coords.90<-Rotation(coords,90*pi/180)
coords.90

plot(coords,xlim=range(rbind(coords.90,coords)[,1]),ylim=range(rbind(coords.90,coords)[,2]),asp=1)
points(coords.90,pch=19)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("sacsarlm")
### * sacsarlm

flush(stderr()); flush(stdout())

### Name: sacsarlm
### Title: Spatial simultaneous autoregressive SAC model estimation
### Aliases: sacsarlm
### Keywords: spatial

### ** Examples

data(oldcol)
COL.sacW.eig <- sacsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, 
 nb2listw(COL.nb, style="W"))
summary(COL.sacW.eig, correlation=TRUE)
W <- as(as_dgRMatrix_listw(nb2listw(COL.nb, style="W")), "CsparseMatrix")
trMatc <- trW(W, type="mult")
summary(impacts(COL.sacW.eig, tr=trMatc, R=2000), zstats=TRUE, short=TRUE)
COL.msacW.eig <- sacsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, 
 nb2listw(COL.nb, style="W"), type="sacmixed")
summary(COL.msacW.eig, correlation=TRUE)
summary(impacts(COL.msacW.eig, tr=trMatc, R=2000), zstats=TRUE, short=TRUE)



cleanEx()
nameEx("set.spChkOption")
### * set.spChkOption

flush(stderr()); flush(stdout())

### Name: set.spChkOption
### Title: Control checking of spatial object IDs
### Aliases: set.spChkOption get.spChkOption chkIDs spNamedVec
###   set.VerboseOption get.VerboseOption set.ClusterOption
###   get.ClusterOption set.ZeroPolicyOption get.ZeroPolicyOption
###   set.rlecuyerSeedOption get.rlecuyerSeedOption
### Keywords: spatial

### ** Examples

data(oldcol)
rownames(COL.OLD)
data(columbus)
rownames(columbus)
get.spChkOption()
oldChk <- set.spChkOption(TRUE)
get.spChkOption()
chkIDs(COL.OLD, nb2listw(COL.nb))
chkIDs(columbus, nb2listw(col.gal.nb))
chkIDs(columbus, nb2listw(COL.nb))
tmp <- try(moran.test(spNamedVec("CRIME", COL.OLD), nb2listw(COL.nb)))
print(tmp)
tmp <- try(moran.test(spNamedVec("CRIME", columbus), nb2listw(col.gal.nb)))
print(tmp)
tmp <- try(moran.test(spNamedVec("CRIME", columbus), nb2listw(COL.nb)))
print(tmp)
set.spChkOption(FALSE)
get.spChkOption()
moran.test(spNamedVec("CRIME", columbus), nb2listw(COL.nb))
tmp <- try(moran.test(spNamedVec("CRIME", columbus), nb2listw(COL.nb),
 spChk=TRUE))
print(tmp)
set.spChkOption(oldChk)
get.spChkOption()



cleanEx()
nameEx("similar.listw")
### * similar.listw

flush(stderr()); flush(stdout())

### Name: similar.listw
### Title: Create symmetric similar weights lists
### Aliases: similar.listw
### Keywords: spatial

### ** Examples

data(oldcol)
COL.W <- nb2listw(COL.nb, style="W")
COL.S <- nb2listw(COL.nb, style="S")
sum(log(1 - 0.5 * eigenw(COL.W)))
sum(log(1 - 0.5 * eigenw(similar.listw(COL.W))))
W_J <- as(as_dsTMatrix_listw(similar.listw(COL.W)), "CsparseMatrix")
I <- as_dsCMatrix_I(dim(W_J)[1])
c(determinant(I - 0.5 * W_J, logarithm=TRUE)$modulus)
sum(log(1 - 0.5 * eigenw(COL.S)))
sum(log(1 - 0.5 * eigenw(similar.listw(COL.S))))
W_J <- as(as_dsTMatrix_listw(similar.listw(COL.S)), "CsparseMatrix")
c(determinant(I - 0.5 * W_J, logarithm=TRUE)$modulus)



cleanEx()
nameEx("skater")
### * skater

flush(stderr()); flush(stdout())

### Name: skater
### Title: Spatial 'K'luster Analysis by Tree Edge Removal
### Aliases: skater
### Keywords: cluster tree

### ** Examples

### loading data
bh <- readShapePoly(system.file("etc/shapes/bhicv.shp",
      package="spdep")[1])
### data standardized 
dpad <- data.frame(scale(bh@data[,5:8]))

### neighboorhod list
bh.nb <- poly2nb(bh)

### calculating costs
lcosts <- nbcosts(bh.nb, dpad)

### making listw
nb.w <- nb2listw(bh.nb, lcosts, style="B")

### find a minimum spanning tree
mst.bh <- mstree(nb.w,5)

### the mstree plot
par(mar=c(0,0,0,0))
plot(mst.bh, coordinates(bh), col=2,       
     cex.lab=.7, cex.circles=0.035, fg="blue")
plot(bh, border=gray(.5), add=TRUE)

### three groups with no restriction
res1 <- skater(mst.bh[,1:2], dpad, 2)

### thee groups with minimum population 
res2 <- skater(mst.bh[,1:2], dpad, 2, 200000, bh@data$Pop)

### thee groups with minimun number of areas
res3 <- skater(mst.bh[,1:2], dpad, 2, 3, rep(1,nrow(bh@data)))

### groups frequency
table(res1$groups)
table(res2$groups)
table(res3$groups)

### the skater plot
par(mar=c(0,0,0,0))
plot(res1, coordinates(bh), cex.circles=0.035, cex.lab=.7)

### more one partition
res1b <- skater(res1, dpad, 1)

### length groups frequency
table(res1$groups)
table(res1b$groups)

### the skater plot, using other colors
plot(res1b, coordinates(bh), cex.circles=0.035, cex.lab=.7,
     groups.colors=colors()[(1:length(res1b$ed))*10])

### the Spatial Polygons plot
plot(bh, col=heat.colors(4)[res1b$groups])




graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("sp.correlogram")
### * sp.correlogram

flush(stderr()); flush(stdout())

### Name: sp.correlogram
### Title: Spatial correlogram
### Aliases: sp.correlogram plot.spcor print.spcor
### Keywords: spatial

### ** Examples

example(nc.sids)
ft.SID74 <- sqrt(1000)*(sqrt(nc.sids$SID74/nc.sids$BIR74) +
  sqrt((nc.sids$SID74+1)/nc.sids$BIR74))
tr.SIDS74 <- ft.SID74*sqrt(nc.sids$BIR74)
cspc <- sp.correlogram(ncCC89_nb, tr.SIDS74, order=8, method="corr",
 zero.policy=TRUE)
print(cspc)
plot(cspc)
Ispc <- sp.correlogram(ncCC89_nb, tr.SIDS74, order=8, method="I",
 zero.policy=TRUE)
print(Ispc)
print(Ispc, "bonferroni")
plot(Ispc)
Cspc <- sp.correlogram(ncCC89_nb, tr.SIDS74, order=8, method="C",
 zero.policy=TRUE)
print(Cspc)
print(Cspc, "bonferroni")
plot(Cspc)
drop.no.neighs <- !(1:length(ncCC89_nb) %in% which(card(ncCC89_nb) == 0))
sub.ncCC89.nb <- subset(ncCC89_nb, drop.no.neighs)
plot(sp.correlogram(sub.ncCC89.nb, subset(tr.SIDS74,  drop.no.neighs),
 order=8, method="corr"))



cleanEx()
nameEx("sp.mantel.mc")
### * sp.mantel.mc

flush(stderr()); flush(stdout())

### Name: sp.mantel.mc
### Title: Mantel-Hubert spatial general cross product statistic
### Aliases: sp.mantel.mc plot.mc.sim
### Keywords: spatial

### ** Examples

data(oldcol)
sim1 <- sp.mantel.mc(COL.OLD$CRIME, nb2listw(COL.nb),
 nsim=99, type="geary", alternative="less")
sim1
plot(sim1)
sp.mantel.mc(COL.OLD$CRIME, nb2listw(COL.nb), nsim=99,
 type="sokal", alternative="less")
sp.mantel.mc(COL.OLD$CRIME, nb2listw(COL.nb), nsim=99,
 type="moran")



cleanEx()
nameEx("spautolm")
### * spautolm

flush(stderr()); flush(stdout())

### Name: spautolm
### Title: Spatial conditional and simultaneous autoregression model
###   estimation
### Aliases: spautolm residuals.spautolm deviance.spautolm coef.spautolm
###   fitted.spautolm print.spautolm summary.spautolm LR1.spautolm
###   logLik.spautolm print.summary.spautolm
### Keywords: spatial

### ** Examples

example(NY_data)
lm0 <- lm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata)
summary(lm0)
lm0w <- lm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata, weights=POP8)
summary(lm0w)
esar0 <- errorsarlm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
 listw=listw_NY)
summary(esar0)
system.time(esar1f <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
 data=nydata, listw=listw_NY, family="SAR", method="full", verbose=TRUE))
summary(esar1f)
system.time(esar1M <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
 data=nydata, listw=listw_NY, family="SAR", method="Matrix", verbose=TRUE))
summary(esar1M)
## Not run: 
##D system.time(esar1M <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
##D  data=nydata, listw=listw_NY, family="SAR", method="Matrix", verbose=TRUE,
##D  control=list(super=TRUE)))
##D summary(esar1M)
##D system.time(esar1s <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
##D  data=nydata, listw=listw_NY, family="SAR", method="spam", verbose=TRUE))
##D summary(esar1s)
##D esar1wf <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
##D  listw=listw_NY, weights=POP8, family="SAR", method="full")
##D summary(esar1wf)
##D system.time(esar1wM <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
##D  data=nydata, listw=listw_NY, weights=POP8, family="SAR", method="Matrix"))
##D summary(esar1wM)
##D esar1ws <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
##D  listw=listw_NY, weights=POP8, family="SAR", method="spam")
##D summary(esar1ws)
##D esar1wlu <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
##D  listw=listw_NY, weights=POP8, family="SAR", method="LU")
##D summary(esar1wlu)
##D esar1wch <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
##D  listw=listw_NY, weights=POP8, family="SAR", method="Chebyshev")
##D summary(esar1wch)
##D ecar1f <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
##D  listw=listw_NY, family="CAR", method="full")
##D summary(ecar1f)
##D system.time(ecar1M <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
##D  data=nydata, listw=listw_NY, family="CAR", method="Matrix"))
##D summary(ecar1M)
##D ecar1s <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
##D  listw=listw_NY, family="CAR", method="spam")
##D summary(ecar1s)
##D ecar1wf <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
##D  listw=listw_NY, weights=nydata$POP8, family="CAR", method="full")
##D summary(ecar1wf)
##D system.time(ecar1wM <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME,
##D  data=nydata, listw=listw_NY, weights=POP8, family="CAR", method="Matrix"))
##D summary(ecar1wM)
##D ecar1ws <- spautolm(Z ~ PEXPOSURE + PCTAGE65P + PCTOWNHOME, data=nydata,
##D  listw=listw_NY, weights=POP8, family="CAR", method="spam")
##D summary(ecar1ws)
##D example(nc.sids)
##D ft.SID74 <- sqrt(1000)*(sqrt(nc.sids$SID74/nc.sids$BIR74) + sqrt((nc.sids$SID74+1)/nc.sids$BIR74))
##D lm_nc <- lm(ft.SID74 ~ 1)
##D sids.nhbr30 <- dnearneigh(cbind(nc.sids$east, nc.sids$north), 0, 30, row.names=row.names(nc.sids))
##D sids.nhbr30.dist <- nbdists(sids.nhbr30, cbind(nc.sids$east, nc.sids$north))
##D sids.nhbr <- listw2sn(nb2listw(sids.nhbr30, glist=sids.nhbr30.dist, style="B", zero.policy=TRUE))
##D dij <- sids.nhbr[,3]
##D n <- nc.sids$BIR74
##D el1 <- min(dij)/dij
##D el2 <- sqrt(n[sids.nhbr$to]/n[sids.nhbr$from])
##D sids.nhbr$weights <- el1*el2
##D sids.nhbr.listw <- sn2listw(sids.nhbr)
##D both <- factor(paste(nc.sids$L_id, nc.sids$M_id, sep=":"))
##D ft.NWBIR74 <- sqrt(1000)*(sqrt(nc.sids$NWBIR74/nc.sids$BIR74) + sqrt((nc.sids$NWBIR74+1)/nc.sids$BIR74))
##D mdata <- data.frame(both, ft.NWBIR74, ft.SID74, BIR74=nc.sids$BIR74)
##D outl <- which.max(rstandard(lm_nc))
##D as.character(nc.sids$names[outl])
##D mdata.4 <- mdata[-outl,]
##D W <- listw2mat(sids.nhbr.listw)
##D W.4 <- W[-outl, -outl]
##D sids.nhbr.listw.4 <- mat2listw(W.4)
##D esarI <- errorsarlm(ft.SID74 ~ 1, data=mdata, listw=sids.nhbr.listw,
##D  zero.policy=TRUE)
##D summary(esarI)
##D esarIa <- spautolm(ft.SID74 ~ 1, data=mdata, listw=sids.nhbr.listw,
##D  family="SAR")
##D summary(esarIa)
##D esarIV <- errorsarlm(ft.SID74 ~ ft.NWBIR74, data=mdata, listw=sids.nhbr.listw,
##D  zero.policy=TRUE)
##D summary(esarIV)
##D esarIVa <- spautolm(ft.SID74 ~ ft.NWBIR74, data=mdata, listw=sids.nhbr.listw,
##D  family="SAR")
##D summary(esarIVa)
##D esarIaw <- spautolm(ft.SID74 ~ 1, data=mdata, listw=sids.nhbr.listw,
##D  weights=BIR74, family="SAR")
##D summary(esarIaw)
##D esarIIaw <- spautolm(ft.SID74 ~ both - 1, data=mdata, listw=sids.nhbr.listw,
##D  weights=BIR74, family="SAR")
##D summary(esarIIaw)
##D esarIVaw <- spautolm(ft.SID74 ~ ft.NWBIR74, data=mdata,
##D  listw=sids.nhbr.listw, weights=BIR74, family="SAR")
##D summary(esarIVaw)
##D ecarIaw <- spautolm(ft.SID74 ~ 1, data=mdata.4, listw=sids.nhbr.listw.4,
##D  weights=BIR74, family="CAR")
##D summary(ecarIaw)
##D ecarIIaw <- spautolm(ft.SID74 ~ both - 1, data=mdata.4,
##D  listw=sids.nhbr.listw.4, weights=BIR74, family="CAR")
##D summary(ecarIIaw)
##D ecarIVaw <- spautolm(ft.SID74 ~ ft.NWBIR74, data=mdata.4,
##D  listw=sids.nhbr.listw.4, weights=BIR74, family="CAR")
##D summary(ecarIVaw)
##D nc.sids$fitIV <- append(fitted.values(ecarIVaw), NA, outl-1)
##D spplot(nc.sids, c("fitIV"), cuts=12) # Cressie 1993, p. 565
##D data(oldcol)
##D COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
##D  nb2listw(COL.nb, style="W"))
##D summary(COL.errW.eig)
##D COL.errW.sar <- spautolm(CRIME ~ INC + HOVAL, data=COL.OLD,
##D  nb2listw(COL.nb, style="W"))
##D summary(COL.errW.sar)
##D data(boston)
##D gp1 <- spautolm(log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2)
##D  + I(RM^2) + AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT), 
##D  data=boston.c, nb2listw(boston.soi), family="SMA")
##D summary(gp1)
## End(Not run)



cleanEx()
nameEx("spweights.constants")
### * spweights.constants

flush(stderr()); flush(stdout())

### Name: spweights.constants
### Title: Provides constants for spatial weights matrices
### Aliases: spweights.constants Szero
### Keywords: spatial

### ** Examples

data(oldcol)
B <- spweights.constants(nb2listw(COL.nb, style="B"))
W <- spweights.constants(nb2listw(COL.nb, style="W"))
C <- spweights.constants(nb2listw(COL.nb, style="C"))
S <- spweights.constants(nb2listw(COL.nb, style="S"))
U <- spweights.constants(nb2listw(COL.nb, style="U"))
print(data.frame(rbind(unlist(B), unlist(W), unlist(C), unlist(S), unlist(U)),
  row.names=c("B", "W", "C", "S", "U")))



cleanEx()
nameEx("ssw")
### * ssw

flush(stderr()); flush(stdout())

### Name: ssw
### Title: Compute the sum of dissimilarity
### Aliases: ssw
### Keywords: cluster multivariate

### ** Examples

data(USArrests)
n <- nrow(USArrests)
ssw(USArrests, 1:n)
ssw(USArrests, 1:(n/2))
ssw(USArrests, (n/2+1):n)
ssw(USArrests, 1:(n/2)) + ssw(USArrests, (n/2+1):n)



cleanEx()
nameEx("stsls")
### * stsls

flush(stderr()); flush(stdout())

### Name: stsls
### Title: Generalized spatial two stage least squares
### Aliases: stsls print.stsls print.summary.stsls summary.stsls
###   residuals.stsls coef.stsls deviance.stsls
### Keywords: spatial

### ** Examples

data(oldcol)
COL.lag.eig <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD, nb2listw(COL.nb))
summary(COL.lag.eig, correlation=TRUE)
COL.lag.stsls <- stsls(CRIME ~ INC + HOVAL, data=COL.OLD, nb2listw(COL.nb))
summary(COL.lag.stsls, correlation=TRUE)
COL.lag.stslsW <- stsls(CRIME ~ INC + HOVAL, data=COL.OLD, nb2listw(COL.nb), W2X=FALSE)
summary(COL.lag.stslsW, correlation=TRUE)
COL.lag.stslsR <- stsls(CRIME ~ INC + HOVAL, data=COL.OLD, nb2listw(COL.nb),
robust=TRUE, W2X=FALSE)
summary(COL.lag.stslsR, correlation=TRUE)
COL.lag.stslsRl <- stsls(CRIME ~ INC + HOVAL, data=COL.OLD, nb2listw(COL.nb),
robust=TRUE, legacy=TRUE, W2X=FALSE)
summary(COL.lag.stslsRl, correlation=TRUE)
data(boston)
gp2a <- stsls(log(CMEDV) ~ CRIM + ZN + INDUS + CHAS + I(NOX^2) + I(RM^2) +
  AGE + log(DIS) + log(RAD) + TAX + PTRATIO + B + log(LSTAT),
 data=boston.c, nb2listw(boston.soi))
summary(gp2a)



cleanEx()
nameEx("subset.listw")
### * subset.listw

flush(stderr()); flush(stdout())

### Name: subset.listw
### Title: Subset a spatial weights list
### Aliases: subset.listw
### Keywords: spatial

### ** Examples

example(columbus)
to.be.dropped <- c(31, 34, 36, 39, 42, 46)
pre <- nb2listw(col.gal.nb)
print(pre)
post <- subset(pre, !(1:length(col.gal.nb) %in% to.be.dropped))
print(post)



cleanEx()
nameEx("subset.nb")
### * subset.nb

flush(stderr()); flush(stdout())

### Name: subset.nb
### Title: Subset a neighbours list
### Aliases: subset.nb
### Keywords: spatial

### ** Examples

example(columbus)
coords <- coordinates(columbus)
plot(col.gal.nb, coords)
to.be.dropped <- c(31, 34, 36, 39, 42, 46)
text(coords[to.be.dropped,1], coords[to.be.dropped,2], labels=to.be.dropped,
  pos=2, offset=0.3)
sub.col.gal.nb <- subset(col.gal.nb,
  !(1:length(col.gal.nb) %in% to.be.dropped))
plot(sub.col.gal.nb, coords[-to.be.dropped,], col="red", add=TRUE)
which(!(attr(col.gal.nb, "region.id") %in%
  attr(sub.col.gal.nb, "region.id")))



cleanEx()
nameEx("summary.nb")
### * summary.nb

flush(stderr()); flush(stdout())

### Name: summary.nb
### Title: Print and summary function for neighbours and weights lists
### Aliases: summary.nb print.nb summary.listw print.listw
### Keywords: spatial

### ** Examples

example(columbus)
coords <- coordinates(columbus)
col.gal.nb
summary(col.gal.nb, coords)
col.listw <- nb2listw(col.gal.nb, style="W")
col.listw
summary(col.listw)



cleanEx()
nameEx("summary.sarlm")
### * summary.sarlm

flush(stderr()); flush(stdout())

### Name: summary.sarlm
### Title: summary method for class sarlm
### Aliases: summary.sarlm print.sarlm print.summary.sarlm
### Keywords: spatial

### ** Examples

data(oldcol)
COL.mix.eig <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb), type="mixed", method="eigen")
summary(COL.mix.eig, correlation=TRUE, Nagelkerke=TRUE)
COL.mix.M <- lagsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
 nb2listw(COL.nb), type="mixed", method="Matrix")
summary(COL.mix.M, correlation=TRUE, Nagelkerke=TRUE)
COL.errW.eig <- errorsarlm(CRIME ~ INC + HOVAL, data=COL.OLD,
  nb2listw(COL.nb, style="W"), method="eigen")
summary(COL.errW.eig, correlation=TRUE, Nagelkerke=TRUE, Hausman=TRUE)



cleanEx()
nameEx("testnb")
### * testnb

flush(stderr()); flush(stdout())

### Name: is.symmetric.nb
### Title: Test a neighbours list for symmetry
### Aliases: is.symmetric.nb sym.attr.nb make.sym.nb old.make.sym.nb
###   is.symmetric.glist
### Keywords: spatial

### ** Examples

example(columbus)
coords <- coordinates(columbus)
ind <- sapply(slot(columbus, "polygons"), function(x) slot(x, "ID"))
print(is.symmetric.nb(col.gal.nb, verbose=TRUE, force=TRUE))
k4 <- knn2nb(knearneigh(coords, k=4), row.names=ind)
k4 <- sym.attr.nb(k4)
print(is.symmetric.nb(k4))
k4.sym <- make.sym.nb(k4)
print(is.symmetric.nb(k4.sym))



cleanEx()
nameEx("tolerance.nb")
### * tolerance.nb

flush(stderr()); flush(stdout())

### Name: tolerance.nb
### Title: Function to construct edges based on a tolerance angle and a
###   maximum distance
### Aliases: tolerance.nb
### Keywords: spatial

### ** Examples

set.seed(1)
ex.data<-cbind(runif(50),rexp(50))

### Construct object of class nb with a tolerance angle of 30 degrees and a maximum distance of 2 m.
nb.ex<-tolerance.nb(ex.data, unit.angle = "degrees", max.dist=1, tolerance = 30)

### Construct object of class nb with a tolerance angle of 30 degrees and a maximum distance of 2 m. The coordinates are rotated at an angle of 45 degrees counterclockwise.
nb.ex2<-tolerance.nb(ex.data, unit.angle = "degrees", max.dist=1, tolerance = 30, rot.angle = 45)

### Construct object of class nb with a tolerance angle of pi/8 radians and a maximum distance of 1.5 m. The coordinates are rotated at an angle of pi/4 radians clockwise.
nb.ex3<-tolerance.nb(ex.data, unit.angle = "radians", max.dist=1.5, tolerance = pi/8,rot.angle = -pi*2/3)


par(mfrow=c(1,3))
plot(nb.ex,ex.data,asp=1)
plot(nb.ex2,ex.data,asp=1)
plot(nb.ex3,ex.data,asp=1)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("trW")
### * trW

flush(stderr()); flush(stdout())

### Name: trW
### Title: Spatial weights matrix powers traces
### Aliases: trW
### Keywords: spatial

### ** Examples

example(columbus)
listw <- nb2listw(col.gal.nb)
W <- as(as_dgRMatrix_listw(listw), "CsparseMatrix")
system.time(trMat <- trW(W, type="mult"))
str(trMat)
set.seed(1100)
system.time(trMC <- trW(W, type="MC"))
str(trMC)
plot(trMat, trMC)
abline(a=0, b=1)
for(i in 3:length(trMC)) segments(trMat[i], trMC[i]-2*attr(trMC, "sd")[i], trMat[i], trMC[i]+2*attr(trMC, "sd")[i])
listwS <- similar.listw(listw)
W <- as(as(as_dgRMatrix_listw(listwS), "CsparseMatrix"), "symmetricMatrix")
system.time(trmom <- trW(W, m=24, type="moments"))
str(trmom)
all.equal(trMat[1:24], trmom, check.attributes=FALSE)
system.time(trMat <- trW(W, m=24, type="mult"))
str(trMat)
all.equal(trMat, trmom, check.attributes=FALSE)
set.seed(1)
system.time(trMC <- trW(W, m=24, type="MC"))
str(trMC)



cleanEx()
nameEx("tri2nb")
### * tri2nb

flush(stderr()); flush(stdout())

### Name: tri2nb
### Title: Neighbours list from tri object
### Aliases: tri2nb
### Keywords: spatial

### ** Examples

example(columbus)
coords <- coordinates(columbus)
ind <- sapply(slot(columbus, "polygons"), function(x) slot(x, "ID"))
col.tri.nb <- tri2nb(coords, row.names=ind)
plot(columbus, border="grey")
plot(col.tri.nb, coords, add=TRUE)
title(main="Raw triangulation links")
x <- seq(0,1,0.1)
y <- seq(0,2,0.2)
xy <- expand.grid(x, y)
try(xy.nb <- tri2nb(xy))
seed <- 1234
xid <- sample(1:nrow(xy))
xy.nb <- tri2nb(xy[xid,])
plot(xy.nb, xy[xid,])



cleanEx()
nameEx("used.cars")
### * used.cars

flush(stderr()); flush(stdout())

### Name: used.cars
### Title: US 1960 used car prices
### Aliases: used.cars usa48.nb
### Keywords: datasets

### ** Examples

data(used.cars)
moran.test(used.cars$price.1960, nb2listw(usa48.nb))
moran.plot(used.cars$price.1960, nb2listw(usa48.nb),
  labels=rownames(used.cars))
uc.lm <- lm(price.1960 ~ tax.charges, data=used.cars)
summary(uc.lm)
lm.morantest(uc.lm, nb2listw(usa48.nb))
lm.morantest.sad(uc.lm, nb2listw(usa48.nb))
lm.LMtests(uc.lm, nb2listw(usa48.nb))
uc.err <- errorsarlm(price.1960 ~ tax.charges, data=used.cars,
  nb2listw(usa48.nb), tol.solve=1.0e-13, control=list(tol.opt=.Machine$double.eps^0.3))
summary(uc.err)
uc.lag <- lagsarlm(price.1960 ~ tax.charges, data=used.cars,
  nb2listw(usa48.nb), tol.solve=1.0e-13, control=list(tol.opt=.Machine$double.eps^0.3))
summary(uc.lag)
uc.lag1 <- lagsarlm(price.1960 ~ 1, data=used.cars,
  nb2listw(usa48.nb), tol.solve=1.0e-13, control=list(tol.opt=.Machine$double.eps^0.3))
summary(uc.lag1)
uc.err1 <- errorsarlm(price.1960 ~ 1, data=used.cars,
  nb2listw(usa48.nb), tol.solve=1.0e-13, control=list(tol.opt=.Machine$double.eps^0.3))
summary(uc.err1)



cleanEx()
nameEx("wheat")
### * wheat

flush(stderr()); flush(stdout())

### Name: wheat
### Title: Mercer and Hall wheat yield data
### Aliases: wheat
### Keywords: datasets

### ** Examples

## Not run: 
##D data(wheat)
##D names(wheat) <- c('lon','lat','yield')
##D wheat$lat1 <- 69 - wheat$lat
##D wheat$r <- factor(wheat$lat1)
##D wheat$c <- factor(wheat$lon)
##D wheat_sp <- wheat
##D coordinates(wheat_sp) <- c("lon", "lat1")
##D wheat_spg <- wheat_sp
##D gridded(wheat_spg) <- TRUE
##D wheat_spl <- as(wheat_spg, "SpatialPolygons")
##D df <- as(wheat_spg, "data.frame")
##D row.names(df) <- sapply(slot(wheat_spl, "polygons"),
##D  function(x) slot(x, "ID"))
##D wheat <- SpatialPolygonsDataFrame(wheat_spl, data=df)
## End(Not run)
wheat <- readShapeSpatial(system.file("etc/shapes/wheat.shp",
 package="spdep")[1])



cleanEx()
nameEx("write.nb.gal")
### * write.nb.gal

flush(stderr()); flush(stdout())

### Name: write.nb.gal
### Title: Write a neighbours list as a GAL lattice file
### Aliases: write.nb.gal
### Keywords: spatial

### ** Examples

example(columbus)
GALfile <- tempfile("GAL")
write.nb.gal(col.gal.nb, GALfile)
col.queen <- read.gal(GALfile)
summary(diffnb(col.queen, col.gal.nb))



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
