aple.mc <- function(x, listw, nsim, override_similarity_check=FALSE) {
    aple.boot <- function(var, i, ...) {
        var <- var[i]
        return(inAple(x=var, ...))
    }
    pre <- preAple(x=x, listw=listw,
        override_similarity_check=override_similarity_check)
    res <- boot(x, statistic=aple.boot, R=nsim, sim="permutation", pre=pre)
    res
}
