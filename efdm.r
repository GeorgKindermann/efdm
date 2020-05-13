require(distances)

require(Matrix)
require(MASS)
require(parallel)


c.factor <- function(..., recursive=TRUE) unlist(list(...), recursive=recursive)

efdmBreakNamesGet <- function(x, env=breaks) {
  i <- do.call(rbind, lapply(ls(breaks), function(y) {
    i <- grep(paste0("^", y), x)
    if(length(i)>0) {data.frame(i, y, stringsAsFactors = FALSE)}}))
  if(length(x) == nrow(i)) {i[order(i[,1]),]}
}

efdmClassGet <- function(x, int=TRUE, env=breaks) {
  i <- efdmBreakNamesGet(colnames(x), env)
  as.data.frame(lapply(seq_len(nrow(i)), function(j) {
    k <- findInterval(x[,j], get(i[j,2], env))
    if(int) {k
    } else {factor(k,0:length(get(i[j,2], env)))}
  }), col.names=colnames(x))
}

efdmDimGet <- function(x, env=breaks) {
  vapply(unique(efdmBreakNamesGet(x, env)[,2]), function(i) {length(get(i, env))}, numeric(1))+1
}

efdmIndexGet <- function(x, d) {
  i <- lapply(paste0("^", names(d)), grep, colnames(x))
  as.vector(matrix(unlist(x), nrow(x)) %*% cumprod(c(1,d))[match(seq_len(ncol(x)), i)]) + 1
}

efdmIndexInv <- function(i, d) {
  m <- unname(rev(cumprod(d))[-1])
  "colnames<-"(do.call(rbind, lapply(i, function(i) rev(Reduce("%%", m, (i-1), accumulate = TRUE) %/% c(m, 1)))), names(d))
}

efdmIndexStatic <- function(i, d) {
  if(length(i) == 0) {0
  } else {
    if(is.character(i)) {i <- match(i, names(d))}
    m <- c(1, cumprod(d))
    rowSums(expand.grid(lapply(i, function(i) seq(0, by=m[i], length=d[i]))))
  }
}

efdmStateGet <- function(x, what="area", breaks=breaks, sparse=TRUE) {
  i <- colnames(x) == what
  t0 <- efdmDimGet(colnames(x)[!i], breaks)
  if(sparse) {   #Sparse array
    x <- aggregate(x[,i]
           , list(i=efdmIndexGet(efdmClassGet(x[!i], TRUE, breaks), t0)), sum)
    x <- Matrix::sparseVector(x$x, x$i, prod(t0))
  } else {     #Dense array
    x <- as.vector(xtabs(reformulate(termlabels = ".", response = what), data.frame(x[i], efdmClassGet(x[!i], FALSE, breaks))))
  }
  attr(x, "dimS") <- t0
  x
}

efdmTransitionGet <- function(x, area="area", flow=c("harvest","residuals"), t0="0", t1="1", breaks=breaks, model=c("qda","lda"), dimWgt=1) {
  env <- new.env(parent=emptyenv())
  env$dimS <- efdmDimGet(colnames(x)[endsWith(colnames(x), t0)], breaks)
  env$dimWgt <- rep_len(dimWgt, length(env$dimS))
  x <- aggregate(. ~ from + to, data.frame(from = efdmIndexGet(efdmClassGet(x[,endsWith(colnames(x), t0)], TRUE, breaks), env$dimS)
   , to = efdmIndexGet(efdmClassGet(x[,endsWith(colnames(x), t1)], TRUE, breaks), env$dimS)
   , tmp[match(c(area,flow), colnames(tmp))]), sum)
  tmp <- aggregate(area~from, x, sum)
  x$share <- x$area / tmp$area[match(x$from, tmp$from)]
  x[flow] <- x[flow] / x[,area]
  env$simple <- x[c("from","to","share",flow)]
  if(!is.null(model) && !is.na(model[1])) {
    fun <- function(i, j) {
      match.fun(model[j])(reformulate(termlabels = names(env$dimS)
                , response = paste0(i,".1")), tt, weights= x[area])}
    tt <- data.frame(efdmIndexInv(x$from,env$dimS),efdmIndexInv(x$to,env$dimS))
    env$model <- lapply(setNames(names(env$dimS), names(env$dimS)), function(i){
      for(j in seq_along(model)) {
        try(ret <- fun(i,j), TRUE)
        if(exists("ret", inherits = FALSE)) return(ret)
      }})
  }
  env
}

efdmNextArea <- function(t0, transition, share="share", nMinToGo=1, nMaxToGo=3, pMinToGo=0.05, lowMem=FALSE, maxNeighbours=7) {
  t1 <- if(is.vector(t0)) {numeric(length(t0))
        } else {Matrix::sparseVector(0, 1, length(t0))}
  attr(t1, "dimS") <- attr(t0, "dimS")
  static <- efdmIndexStatic(setdiff(names(attr(t0, "dimS")), names(transition$dimS)), attr(t0, "dimS"))
  from <- efdmIndexGet(efdmIndexInv(transition$simple$from, transition$dimS), attr(t0, "dimS"))
  to <- efdmIndexGet(efdmIndexInv(transition$simple$to, transition$dimS), attr(t0, "dimS"))
  tt <- aggregate(as.vector(t0[rowSums(expand.grid(from, static))]) * transition$simple[,share], list(to=rowSums(expand.grid(to, static))), sum)
  tt <- tt[tt$x > 0,]
  t1[tt$to] <- as.vector(t1[tt$to]) + tt$x
  i <- if(is.vector(t0)) {which(t0 != 0)} else {t0@i}
  i <- i[!i %in% rowSums(expand.grid(from, static))]
  if(length(i) > 0) {
    if(exists("model", envir = transition)) {
      if(lowMem) {
        lapply(i, function(i) {
          from <- as.data.frame(efdmIndexInv(i, attr(t0, "dimS")))
          to <- lapply(setNames(seq_along(transition$model), names(transition$model)), function(i) {
            f <- transition$model[[i]]
            if(is.null(f)) {setNames(1,from[,names(transition$model)[i]])
            } else {
              x <- predict(f, from)$posterior
              names(x) <- colnames(x)
              proportions(head(x[order(-x)], min(nMaxToGo, max(nMinToGo, sum(x >= pMinToGo)))))}
          })
          toi <- efdmIndexGet(expand.grid(lapply(to, function(x) as.integer(names(x)))), attr(t1, "dimS"))
          t1[toi] <<- as.vector(t1[toi]) + as.vector(t0[i]) * apply(expand.grid(to), 1, prod)
        })
      } else {
        from <- as.data.frame(efdmIndexInv(i, attr(t0, "dimS")))
        to <- lapply(setNames(seq_along(transition$model), names(transition$model)), function(i) {
          f <- transition$model[[i]]
          if(is.null(f)) {setNames(rep(1,nrow(from)),from[,names(transition$model)[i]])
          } else {
            apply(predict(f, from)$posterior, 1, function(x) proportions(head(x[order(-x)], min(nMaxToGo, max(nMinToGo, sum(x >= pMinToGo))))))
          }})
        rm(from)
        x <- do.call(cbind, to)
        to <- do.call(cbind, lapply(to, function(x) if(is.list(x)) lapply(x, function(x) as.integer(names(x))) else as.integer(names(x))))
        x <- apply(x, 1, function(x) apply(expand.grid(x), 1, prod))
        to <- apply(to, 1, function(x) efdmIndexGet(expand.grid(x), attr(t1, "dimS")))
        x <- mapply("*", x, as.vector(t0[i]))
        to <- aggregate(unlist(x), list(to=unlist(to)), sum)
        rm(x)
        t1[to$to] <- as.vector(t1[to$to]) + to$x
      }
    } else if(exists("similar", envir = transition)) {
      from <- as.matrix(as.data.frame(efdmIndexInv(i, attr(t0, "dimS"))))
      from <- from[,names(transition$dimS)]
      fromObsI <- unique(transition$simple$from)
      fromObs <- as.matrix(as.data.frame(efdmIndexInv(fromObsI, transition$dimS)))
      x <- rbind(fromObs, from)
      x <- distances(x * transition$dimWgt[col(x)])
      idx <- seq_len(nrow(from))+nrow(fromObs)
      to <- nearest_neighbor_search(x, maxNeighbours, idx, seq_len(length(fromObsI)))
      to <- lapply(seq_along(idx), function(j) {
        y <- x[idx[j],to[,j]]
        to[,j][y == y[1]]
      })
      rm(x, idx)
#      to <- mclapply(split(from, seq_len(nrow(from))), function(y) {
#        x <- abs(fromObs - y[col(fromObs)]) %*% transition$dimWgt
#        which(x == min(x))
#      })
      transition$simple <- transition$simple[order(transition$simple$from),]
      x <- cbind(1L+findInterval(fromObsI[unlist(to)]-1, transition$simple$from), findInterval(fromObsI[unlist(to)], transition$simple$from))
      x <- mapply(seq, x[,1], x[,2])
      to <- unname(split(x, rep(seq_along(to), lengths(to))))
      rm(x)
      to <- lapply(to, function(y) {j <- unlist(y)
        list(transition$simple$to[j], proportions(transition$simple$share[j])
           , transition$simple$from[j]) #Only needed for relativ shift
      })
#      to <- mclapply(to, function(y) {
#        j <- which(transition$simple$from %in% fromObsI[y])
#        list(transition$simple$to[j], proportions(transition$simple$share[j])
#           , transition$simple$from[j]) #Only needed for relativ shift
#      })
      to <- do.call(rbind, to)
      to[,2] <- mapply("*", to[,2], as.vector(t0[i]))
      if(exists("absolute", envir = transition) & transition$absolute) {
        to[,1] <- mapply(function(x,y) {
          tt <- efdmIndexInv(x, transition$dimS)
          tt[!transition$similar[col(tt)]] <- rep(y[!transition$similar], each=nrow(tt))
          efdmIndexGet(tt, attr(t1, "dimS"))
        }, to[,1], split(from, seq_len(nrow(from))))
      } else {
        to[,1] <- mapply(function(x,y) {
          tt <- x + y[col(x)]
          j <- tt >= transition$dimS[col(tt)]
          tt[j] <- (transition$dimS - 1)[col(tt)][j]
          j <- tt < 0
          tt[j] <- 0
          efdmIndexGet(tt, attr(t1, "dimS"))}, 
          mapply(function(x,y) {
            efdmIndexInv(x, transition$dimS) - efdmIndexInv(y, transition$dimS)},
            to[,1], to[,3]), split(from, seq_len(nrow(from))))
      }
      to <- aggregate(unlist(to[,2]), list(to=unlist(to[,1])), sum)
      t1[to$to] <- as.vector(t1[to$to]) + to$x
    } else if(exists("shifter", envir = transition)) {
      from <- as.matrix(as.data.frame(efdmIndexInv(i, attr(t0, "dimS"))))
      from <- from[,names(transition$dimS)]
      to <- from + transition$shifter[col(from)]
      j <- to >= transition$dimS[col(to)]
      to[j] <- (transition$dimS - 1)[col(to)][j]
      j <- to < 0
      to[j] <- 0
      to <- efdmIndexGet(to, attr(t1, "dimS"))
      to <- aggregate(as.vector(t0[i]), list(to=to), sum)
      t1[to$to] <- as.vector(t1[to$to]) + to$x
    } else {t1[i] <- as.vector(t1[i]) + t0[i]}
  }
  t1
}

efdmNextFlow <- function(t0, transition, flow="harvest", share="share") {
  res <- if(is.vector(t0)) {numeric(length(t0))
         } else {Matrix::sparseVector(0, 1, length(t0))}
  attr(res, "dimS") <- attr(t0, "dimS")
  static <- efdmIndexStatic(setdiff(names(attr(t0, "dimS")), names(transition$dimS)), attr(t0, "dimS"))
  from <- efdmIndexGet(efdmIndexInv(transition$simple$from, transition$dimS), attr(t0, "dimS"))
  tt <- aggregate(as.vector(t0[rowSums(expand.grid(from, static))]) * transition$simple[,share] * transition$simple[,flow], list(to=rowSums(expand.grid(from, static))), sum)
  tt <- tt[tt$x > 0,]
  res[tt$to] <- tt$x
  i <- if(is.vector(t0)) {which(t0 != 0)} else {t0@i}
  i <- i[!i %in% rowSums(expand.grid(from, static))]
  if(length(i) > 0) {
    tt <- aggregate(cbind(x = transition$simple[,flow] * transition$simple[,share]), list(from=transition$simple$from), sum)
    tFrom <- efdmIndexInv(tt$from, transition$dimS)
    iFrom <- efdmIndexInv(i, attr(t0, "dimS"))
    iFrom <- iFrom[,match(colnames(tFrom), colnames(iFrom))]
    res[i] <- res[i] + t0[i] * vapply(seq_along(i), function(k) {
      j <- abs(sweep(tFrom, 2, iFrom[k,])) %*% transition$dimWgt
      mean(tt$x[which(j == min(j))])}, numeric(1))
  }
  res
}
