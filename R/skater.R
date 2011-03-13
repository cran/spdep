skater <- function(edges, data, ncuts, crit, vec.crit,
                   method=c("euclidean", "maximum", "manhattan",
                     "canberra", "binary", "minkowski",
                     "mahalanobis", "other"), p=2, 
                   cov, inverted=FALSE, otherfun) {
  n <- nrow(data)
  if (!missing(ncuts))
    if (ncuts<1 & n<ncuts)
      stop(paste("Number of cuts must be from 1 to", n-1))
  if (any(class(edges)=="skater"))
    res <- edges
  else {
    res <- list(groups=rep(1, n),
                edges.groups=list(list(node=1:n, edge=edges)),
                not.prune=NULL, candidates=1,
                ssto=ssw(data, 1:n, method, p, cov, inverted, otherfun))
    res$ssw <- res$edges.groups[[1]]$ssw <- res$ssto
    tmp <- sort(prunecost(res$edges.groups[[1]]$edge[,1:2, drop=FALSE], 
                          data, method, p, cov, inverted, otherfun), 
                decreasing=TRUE, method='quick', index.return=TRUE)
    res$edges.groups[[1]]$edge =
      cbind(res$edges.groups[[1]]$edge[tmp$ix, ], tmp$x)
    if (missing(crit))
      res$crit <- 1
    else
      res$crit <- crit
    if (missing(vec.crit))
      res$vec.crit <- rep(1,n)
    else
      res$vec.crit <- vec.crit
  }
  cuts <- length(res$edges.groups)
  if (missing(ncuts))
    ncuts <- n-cuts
  else
    ncuts <- ncuts+cuts-1
  if (is.null(res$vec.crit))
    res$vec.crit <- rep(1, n)
  if (is.null(res$crit))
    res$crit <- 1
  npr <- which(sapply(res$edges.groups, function(x)
                      sum(res$vec.crit[x$node])<2*res$crit))
  res$not.prune <- unique(c(res$not.prune, npr))
  res$candidates <- setdiff(1:length(res$edges.groups), res$not.prune)
  repeat {
    if (cuts>ncuts)
      break
    if (length(res$candidates)==0L)
      break
    l.costs.ord <- lapply(res$edges.groups[res$candidates],
                          function(x) x$edge[,3])
    dc <- cbind(id=rep(res$candidates, sapply(l.costs.ord, length)),
                cost=unlist(l.costs.ord),
                idi=unlist(lapply(l.costs.ord, function(x) 1:length(x))))
    dc <- dc[sort(dc[,2], method="quick", decreasing=TRUE,
                  index.return=TRUE)$ix,, drop=FALSE]
    k <- 1
    repeat {
      toprun <- rbind(res$edges.groups[[dc[k,1]]]$edge[dc[k,3],1:2],
                      res$edges.groups[[dc[k,1]]]$edge[-dc[k,3],1:2])
      g.pruned <- prunemst(toprun, only.nodes=FALSE)
      if (any(sapply(g.pruned, function(x)
                     sum(res$vec.crit[x$node]))<res$crit)) {
        id.not <- !is.element(res$candidates, unique(dc[-(1:k),1]))
        res$not.prune <- unique(c(res$not.prune, res$candidates[id.not]))
        res$candidates <- setdiff(1:length(res$edges.groups),
                                  res$not.prune)
        k <- k + 1
        if (k>nrow(dc)) {
          break
        }
      }
      else {
        gc.pruned <- lapply(g.pruned, function(e) {
          if (nrow(e$edge)==0)
            return(list(node=e$node, edge=matrix(0,0,3),
                        ssw=ssw(data, e$node, method, p, cov,
                          inverted, otherfun)))
          else {
            tmp <- sort(prunecost(e$edge[, 1:2, drop=FALSE], data,
                                  method, p, cov, inverted, otherfun),
                        decreasing=TRUE, method='quick', index.return=TRUE) 
            list(node=e$node,
                 edge=cbind(e$edge[tmp$ix, , drop=FALSE], tmp$x),
                 ssw=ssw(data, e$node, method, p, cov, inverted, otherfun))
          }
        })
        res$edges.groups[[dc[k,1]]] <- gc.pruned[[1]]
        cuts <- cuts + 1
        res$edges.groups[[cuts]] <- gc.pruned[[2]]
        res$ssw <- c(res$ssw, sum(sapply(res$edges.groups,
                                         function(e) sum(e$ssw))))
        npr <- which(sapply(res$edges.groups, function(x)
                            sum(res$vec.crit[x$node])<2*res$crit))
        res$not.prune <- unique(c(res$not.prune, npr))
        res$candidates <- setdiff(1:length(res$edges.groups), res$not.prune)
        break
      }
    }
  }
  for (i in 1:length(res$edges.groups))
    res$groups[res$edges.groups[[i]]$node] <- i
  attr(res, "class") <- "skater"
  return(res)
}
