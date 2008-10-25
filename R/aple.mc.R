aple.mc <- function(x, listw, nsim) {
    aple.boot <- function(var, i, ...) {
        var <- var[i]
        return(inAple(x=var, ...))
    }
    pre <- preAple(x=x, listw=listw)
    res <- boot(x, statistic=aple.boot, R=nsim, sim="permutation", pre=pre)
    res
}
