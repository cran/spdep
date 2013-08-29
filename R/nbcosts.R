nbcosts <- function(nb, data, method=c("euclidean", "maximum", "manhattan",
                                "canberra", "binary", "minkowski",
                                "mahalanobis"), p=2, cov, inverted=FALSE) {
  require(parallel)
  clist <- mclapply(1:length(nb), function(i)
                    nbcost(data, i, nb[[i]], method,
                           p, cov, inverted),
                    mc.cores=ifelse(is.null(getOption('mc.cores')),
                        detectCores(), options('mc.cores')))
  attr(clist, "call") <- match.call()
  attr(clist, "class") <- "nbdist"
  return(clist)
}

nbcost <- function(data, id, id.neigh,
                   method=c("euclidean", "maximum", "manhattan",
                     "canberra", "binary", "minkowski",
                     "mahalanobis"), p=2, cov, inverted=FALSE) {
  if (is.function(method))
    return(method(data, id, id.neigh))
  else {
    method <- match.arg(method)
    data <- as.matrix(data)
    if (method=="mahalanobis")
      return(mahalanobis(data[id.neigh,,drop=FALSE], data[id,,drop=FALSE],
cov, inverted))
    else
      return(dist(rbind(data[id,,drop=FALSE], data[id.neigh,,drop=FALSE]),
method=method,
                p=p)[1:length(id.neigh)])
  }
}


