nbcosts <- function(nb, data, method=c("euclidean", "maximum", "manhattan",
                                "canberra", "binary", "minkowski",
                                "mahalanobis", "other"), p=2, 
                    cov, inverted=FALSE, otherfun) {
  clist <- lapply(1:length(nb), function(i)
                  nbcost(data, i, nb[[i]], method,
                         p, cov, inverted, otherfun))
  attr(clist, "call") <- match.call()
  attr(clist, "class") <- "nbdist"
  return(clist)
}

nbcost <- function(data, id, id.neigh,
                   method=c("euclidean", "maximum", "manhattan",
                     "canberra", "binary", "minkowski",
                     "mahalanobis", "other"), p=2, 
                   cov, inverted=FALSE, otherfun) {
  method <- match.arg(method)
  if (method=="other")
    return(otherfun(data, id, id.neigh))
  if (method=="mahalanobis")
    return(mahalanobis(data[id.neigh,], data[id,], cov, inverted))
  else
    return(dist(rbind(data[id,], data[id.neigh, ]), method=method,
                p=p)[1:length(id.neigh)])
}
