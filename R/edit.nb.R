# Copyright 2001-8 by Roger Bivand and Nicholas Lewin-Koh
#

edit.nb <- function(name, coords, polys=NULL, ...) {
  nb <- name
  cnb <- card(nb)
# class to inherits Jari Oksanen 080603
  if (!inherits(x, "nb")) stop("not a neighbours list")
  cl <- class(nb)
  if (length(cl) > 1) icl <- cl[-match("nb", cl)]
  else icl <- NULL
  x <- coords[,1]
  y <- coords[,2]
  n <- length(nb)
  row.names <- attr(nb, "region.id")
  if (is.null(row.names)) row.names <- as.character(1:n)
  xlim <- range(x)
  ylim <- range(y)
  plot.new()
  plot.window(xlim = xlim, ylim = ylim, "", asp=1)
  if (!is.null(polys))
    plot(polys, border="grey", add=TRUE)
  for (i in 1:n) {
    #arrows(x[i],y[i],x[nb[[i]]],y[nb[[i]]],lenght=.08, angle=.15)
###
    if (cnb[i] > 0) segments(x[i],y[i],x[nb[[i]]],y[nb[[i]]])
###
    #inb <- nb[[i]]
    #for (j in inb)
    #lines(c(x[i], x[j]), c(y[i], y[j]), col="black")
  }
  points(x, y)
  finished <- "n"
  deletions <- NULL
  additions <- NULL
###
  edit.segs<-list()
  e.seg.stat<-NULL
  enum<-0
  erase.col<-par()$bg
###
  while (finished == "n") {
    cat("Identifying contiguity for deletion ...\n")
    cand <- identify(x, y, n=2)
    lines(x[cand], y[cand], col="red")

    if (.Platform$OS.type == "windows") bringToTop(-1)
    if ((cand[2] %in% nb[[cand[1]]]) && (cand[1] %in% nb[[cand[2]]])) {
      delete <- readline("Delete this line (y/n) ")
      if (delete != "y") delete <- "n"
      else {
        deletions <- c(deletions, paste(cand, collapse="-"))
        nb[[cand[1]]] <- nb[[cand[1]]][nb[[cand[1]]] != cand[2]]
				if(length(nb[[cand[1]]]) == 0) {
                                  nb[[cand[1]]] <- as.integer(0)
                                  cat(cand[1], "is now an island\n")
				}
        nb[[cand[2]]] <- nb[[cand[2]]][nb[[cand[2]]] != cand[1]]
        if(length(nb[[cand[2]]]) == 0) {
          nb[[cand[2]]] <- as.integer(0)
          cat(cand[2], "is now an island\n")
        }
###
        lines(x[cand], y[cand], col=erase.col)
        lines(x[cand], y[cand], col='brown',lty=4)
        enum<-enum+1
        edit.segs[[enum]]<-cand
        e.seg.stat<-c(e.seg.stat,0)
###
        cat("deleted contiguity between point", cand[1], "and", cand[2], "\n")
      }

      #plot.new()
      #plot.window(xlim = xlim, ylim = ylim, "", asp=1)
      #if (!is.null(polys))
      #  plot(polys, border="grey", add=TRUE)
      #for (i in 1:n) {
      #  inb <- nb[[i]]
      #  for (j in inb)
      #    lines(c(x[i], x[j]), c(y[i], y[j]),
      #          col="black")
      #}
      #points(x, y)
    }
      else {
        if (length(cand == 2)) {
          cat("No contiguity between chosen points\n")
          addcont <- readline("Add contiguity? (y/n) ")
          if (addcont != "y") addcont <- "n"
          if (addcont == "y") {
            nb[[cand[1]]] <-
              sort(unique(c(nb[[cand[1]]], cand[2])))
            nb[[cand[2]]] <-
              sort(unique(c(nb[[cand[2]]], cand[1])))
            cat("added contiguity between point",
                cand[1], "and", cand[2], "\n")
            additions <- c(additions, paste(cand, collapse="-"))
###
            enum<-enum+1
            edit.segs[[enum]]<-cand
            e.seg.stat<-c(e.seg.stat,1)
            lines(x[cand], y[cand], col='yellow')
###
          }
#          plot.new()
#          plot.window(xlim = xlim, ylim = ylim, "", asp=1)
#          if (!is.null(polys))
#            plot(polys, border="grey", add=TRUE)
#          for (i in 1:n) {
#            inb <- nb[[i]]
#            for (j in inb)
#              lines(c(x[i], x[j]),
#                    c(y[i], y[j]),
#                    col="black")
#				}
#          points(x, y)
        }
      }
#    finished <- readline("Finished yet? (y/n) ")
### 
    finished <- readline("Options: quit[q] refresh[r] continue[c] ")
    if (finished == "r") {
      plot.new()
      plot.window(xlim = xlim, ylim = ylim, "", asp=1)
      if (!is.null(polys))
        plot(polys, border="grey", add=TRUE)
      for (i in 1:n) {
        if(nb[[i]][1]!=0 & length(nb[[i]])>0)
          segments(x[i],y[i],x[nb[[i]]],y[nb[[i]]])
      }
      if(enum>1){
        for(i in 1:enum){
          if(e.seg.stat[i]==0){
            lines(x[edit.segs[[i]]], y[edit.segs[[i]]], col=erase.col)
            lines(x[edit.segs[[i]]], y[edit.segs[[i]]], col='brown',lty=4)
          }
          else lines(x[edit.segs[[i]]], y[edit.segs[[i]]], col='yellow')
        }
      }
      points(x, y)
      finished <- readline("Options: quit[q] continue[c]")
    }
    if (finished != "q") finished <- "n"
####
  }
  
  attributes(nb) <- list(deleted=deletions)
  attr(nb, "added") <- additions
  attr(nb, "region.id") <- row.names
  if (is.null(icl)) class(nb) <- "nb"
  else class(nb) <- c("nb", icl)
  nb <- sym.attr.nb(nb)
  nb
}

