prunecost <- function(edges, data,
                      method=c("euclidean", "maximum", "manhattan",
                        "canberra", "binary", "minkowski",
                        "mahalanobis"), p=2, cov, inverted=FALSE) {
  require(parallel)
  sswt <- ssw(data, unique(as.integer(edges)),
              method, p, cov, inverted)
  sswp <- simplify2array(mclapply(1:nrow(edges), function(i) {
    pruned.ids <- prunemst(rbind(edges[i, ], edges[-i, ]),
                           only.nodes=TRUE)
    sum(sapply(pruned.ids, function(j) 
               ssw(data, j, method, p, cov, inverted)))
  }, mc.cores=ifelse(is.null(getOption('mc.cores')),
         detectCores(), options('mc.cores'))))
  return(sswt - sswp)
}

ssw <- function(data, id, method=c("euclidean", "maximum", "manhattan",
                            "canberra", "binary", "minkowski",
                            "mahalanobis"), p=2, cov, inverted=FALSE) {
  if (is.function(method))
    return(method(data, id))
  else { 
    method <- match.arg(method)
    data <- as.matrix(data)
    if (method=="mahalanobis")
      return(sum(mahalanobis(data[id,,drop=FALSE],
                             colMeans(data[id,,drop=FALSE]),
                             cov, inverted)))
    else
      return(sum(dist(rbind(colMeans(data[id,,drop=FALSE]),
                            data[id,,drop=FALSE]),
                      method, p=p)[1:length(id)]))
  }
}

