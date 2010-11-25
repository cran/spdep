aple.mc <- function(x, listw, nsim, override_similarity_check=FALSE,
    useTrace=TRUE) {
    aple.boot <- function(var, i, ...) {
        var <- var[i]
        return(inAple(x=var, ...))
    }
    pre <- preAple(x=x, listw=listw,
        override_similarity_check=override_similarity_check, useTrace=useTrace)
    cl <- get("cl", env = .spdepOptions)
    if (!is.null(cl) && length(cl) > 1) {
        nnsim <- boot_wrapper_in(cl, nsim)
        lres <- clusterCall(cl, boot, x, statistic=aple.boot, R=nnsim,
            sim="permutation", pre=pre)
        res <- boot_wrapper_out(lres, match.call())
    } else {
        res <- boot(x, statistic=aple.boot, R=nsim, sim="permutation", pre=pre)
    }
    res
}

boot_wrapper_in <- function(cl, nsim) {
        require(snow)
        require(rlecuyer)
        rlseed <- get("rlecuyerSeed", env = .spdepOptions)
        if (storage.mode(rlseed) != "integer") rlseed <- as.integer(rlseed)
        if (length(rlseed) != 6) rlseed <- rep(12345, 6)
        clusterSetupRNGstream(cl, seed=rlseed)
        clusterEvalQ(cl, library(boot))
        nnsim <- ceiling(nsim/length(cl))
        nnsim
}

boot_wrapper_out <- function(lres, mcall) {
        res <- list()
        res$t0 <- lres[[1]]$t0
        res$t <- matrix(c(sapply(lres, function(x) x$t)), ncol=1)
        res$R <- sum(sapply(lres, function(x) x$R))
        res$data <- lres[[1]]$data
        res$seed <- c(sapply(lres, function(x) x$seed))
        res$statistic <- lres[[1]]$statistic
        res$sim <- lres[[1]]$sim
        res$call <- mcall
        res$stype <- lres[[1]]$stype
        res$strata <- lres[[1]]$strata
        class(res) <- "boot"
        res
}
